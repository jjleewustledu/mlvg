classdef Test_Fung2013 < matlab.unittest.TestCase
	%% TEST_FUNG2013 

	%  Usage:  >> results = run(mlvg_unittest.Test_Fung2013)
 	%          >> result  = run(mlvg_unittest.Test_Fung2013, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 25-Mar-2021 21:45:21 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlvg/test/+mlvg_unittest.
 	%% It was developed on Matlab 9.9.0.1592791 (R2020b) Update 5 for MACI64.  Copyright 2021 John Joowon Lee.
 	
	properties
        anatPath
        bidsObj
        cacheMat % assign after segmentations, centerlines, registration targets are ready
        corners
        ho
        ho_sumt
        petPath
 		registry
        sourceAnatPath
        sourcePetPath
        t1w
 		testObj
 	end

	methods (Test)
        function test_ctor(this)
            disp(this.testObj)
        end
        function test_cache(this)
            disp(this.testObj)
        end
        function test_t1w(this)
            disp(this.testObj.anatomy)
            this.testObj.anatomy.view()
        end
        function test_buildSegmentation(this)
            this.testObj.buildSegmentation('iterations', 100, 'smoothFactor', 0);
            disp(this.testObj)
        end
		function test_buildCenterlines(this)
            f = this.testObj;        
            f.buildSegmentation('iterations', 130, 'smoothFactor', 0);
            f.buildCenterlines()
        end
        function test_registerCenterline_cpd(this) 
            assert(~isempty(this.cacheMat), 'testObj needs some aufbau to proceed with testing')
            f = this.testObj;          
            f.registerCenterline(f.centerlines_pcs{1}, 'alg', 'cpd', 'laterality', 'L')
            disp(f)
            disp(f.registration)
        end
        function test_registerCenterline_fung(this) 
            assert(~isempty(this.cacheMat), 'testObj needs some aufbau to proceed with testing')
            testObj = this.testObj;          
            testObj.registerCenterline(testObj.centerlines_pcs{1}, 'alg', 'fung', 'laterality', 'L');
            disp(testObj)
            disp(testObj.registration)
        end
        function test_pointCloudToIC(this)
            f = this.testObj;
            f.registerCenterline(f.centerlines_pcs{1}, 'alg', 'cpd', 'laterality', 'L')
            ic = f.pointCloudToIC(f.registration.centerlineOnTarget{1}, 'centerlineOnTarget');
            ic = ic.imdilate(strel('sphere', 2));
            ic.fsleyes
        end
        function test_decay_uncorrected(this)
            obj = this.testObj;
            mask = mlfourd.ImagingContext2('~jjlee/Singularity/CCIR_01211/derivatives/sub-108293/mri/wmparc_on_T1w.nii.gz');
            ho_ = this.ho.volumeAveraged(mask);
            ho_row = obj.decay_uncorrected(ho_);
            plot(obj.timesMid('HO'), ho_row, obj.timesMid('HO'), ho_.nifti.img)
            legend('decay uncorrected', 'decay corrected')
        end        
        function test_call(this)
            [~,ics] = this.testObj.call();
            disp(ics)
%             for i = 1:length(ics)
%                 ics{i}.fsleyes(this.t1w.fqfilename, ics{i}.fqfilename)
%             end
        end
	end

 	methods (TestClassSetup)
		function setupFung2013(this)
 			import mlvg.*;
            this.anatPath = fullfile(getenv('HOME'), 'Singularity/CCIR_01211/derivatives/sub-108293/ses-20210218/anat');
            this.petPath = fullfile(getenv('HOME'), 'Singularity/CCIR_01211/derivatives/sub-108293/ses-20210421/pet');
            this.sourceAnatPath = fullfile(getenv('HOME'), 'Singularity/CCIR_01211/sourcedata/sub-108293/ses-20210218/anat');
            this.sourcePetPath = fullfile(getenv('HOME'), 'Singularity/CCIR_01211/sourcedata/sub-108293/ses-20210421/pet');
            %this.sourcePetPath = fullfile(getenv('HOME'), 'Singularity/subjects/sub-S58163/resampling_restricted');
            cd(this.petPath)         
            %this.corners = [158 122 85; 96 126 88; 156 116 27; 101 113 28]; % PPG
            %this.corners = [113 178 140; 87 178 140; 136 149 58; 62 148 59] + 1; % long vglab
            %this.corners = [140 144 109; 60 144 105; 136 149 58; 62 148 59] + 1; % short vglab; on FreeSurfer T1
            this.corners = [57 149 114; 150 151 114; 66 144 55; 149 145 55] + 1; % T1w vNAV; [ [RS]; [LS]; [RI]; [LI] ].

            this.t1w = mlfourd.ImagingContext2(fullfile(this.sourceAnatPath, 'sub-108293_ses-20210218081506_T1w_MPR_vNav_4e_RMS_orient-std.nii.gz'));
            this.ho = mlfourd.ImagingContext2(fullfile(this.petPath, 'sub-108293_ses-20210421152358_trc-ho_proc-dyn_pet_on_T1w.nii.gz'));
            this.ho_sumt = mlfourd.ImagingContext2(fullfile(this.petPath, 'sub-108293_ses-20210421152358_trc-ho_proc-dyn_pet_on_T1w_avgt.nii.gz'));
            %this.ho = mlfourd.ImagingContext2(fullfile(this.petPath, 'hodt20190523120249_on_T1001.nii.gz'));
            %this.ho_sumt = mlfourd.ImagingContext2(fullfile(this.petPath, 'hodt20190523120249_on_T1001_avgt.nii.gz'));
            %this.cacheMat = fullfile(getenv('HOME'), 'MATLAB-Drive', 'mlvg', 'test', '+mlvg_unittest', 'Test_Fung2013_Vision_20211109.mat');
            %this.cacheMat = fullfile(getenv('HOME'), 'MATLAB-Drive', 'mlvg', 'test', '+mlvg_unittest', 'Test_Fung2013_PPG_20211109.mat');

            this.bidsObj = mlvg.Ccir1211Mediator(this.ho);
 			this.testObj_ = Fung2013(coords=this.corners, plotmore=true, iterations=100, bidsObj=this.bidsObj);
 			%this.testObj_ = Fung2013('coords', this.corners, 'plotmore', true, 'iterations', 1000, 'smoothFactor', 0.1, 'BBBuf', [10 10 4]);




%             b.destinationPath = this.petPath;
%             b.anatPath = ;
%             b.petPath = ;
%             b.t1w_ic = ;
%             b.wmparc_ic = ;
%             b.taus = containers.Map;
%             b.tracer = 'FDG';

 		end
	end

 	methods (TestMethodSetup)
		function setupFung2013Test(this)
 			this.addTeardown(@this.cleanTestMethod);
            if isempty(this.cacheMat)
                this.testObj = copy(this.testObj_);
                return
            end
            if isfile(this.cacheMat)
                ld = load(this.cacheMat, 'testObj');
                this.testObj = ld.testObj;
                return
            end
            this.testObj = copy(this.testObj_);
            this.testObj.buildSegmentation();
            this.testObj.buildCenterlines()                
            this.testObj.buildRegistrationTargets(this.ho)
            if ~isfile(this.cacheMat)
                testObj = this.testObj; %#ok<PROP> 
                save(this.cacheMat, 'testObj')
            end
 		end
	end

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
		function cleanTestMethod(this)
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

