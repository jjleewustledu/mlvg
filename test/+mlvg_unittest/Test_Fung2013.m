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
        cacheMat
        corners
        ho
        ho_sumt
        petPath
 		registry
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
		function test_buildCenterlines(this)
            f = this.testObj;
            f.buildCorners(this.corners);            
            f.buildSegmentation(130, 'smoothFactor', 0);
            f.buildCenterlines()
        end
        function test_registerCenterline(this) 
            f = this.testObj;          
            f.registerCenterline(f.centerlines_pcs{1}, 'alg', 'cpd', 'laterality', 'L')
            disp(f)
            disp(f.registration)            
        end
        function test_pointCloudToIC(this)
            f = this.testObj;
            f.registerCenterline(f.centerlines_pcs{1}, 'alg', 'cpd', 'laterality', 'L')
            ic = f.pointCloudToIC(f.registration.centerlineOnTarget{1}, 'centerlineOnTarget');
            ic = ic.imdilate(strel('sphere', 2));
            ic.fsleyes
        end
        function test_call0(this)            
            this.testObj.buildCorners(this.corners);
            this.testObj.buildSegmentation(100, 'smoothFactor', 0);
            this.testObj.buildCenterlines()
            disp(this.testObj)
        end
        function test_call(this)
            [ics,icrefs] = this.testObj.call();
            disp(ics)
            for i = 1:length(ics)
                ics{i}.fsleyes(this.t1w.fqfilename, icrefs{i}.fqfilename)
            end
        end
	end

 	methods (TestClassSetup)
		function setupFung2013(this)
 			import mlvg.*;
            this.anatPath = '/Users/jjlee/Singularity/CCIR_01211/sourcedata/sub-108293/anat';
            this.petPath = '/Users/jjlee/Singularity/CCIR_01211/sourcedata/sub-108293/pet';
            cd(this.petPath)
            %cd('/Users/jjlee/Singularity/subjects/sub-S58163/resampling_restricted')            
            %this.corners = [113 178 140; 87 178 140; 136 149 58; 62 148 59] + 1; % long
            this.corners = [140 144 109; 60 144 105; 136 149 58; 62 148 59] + 1; % short
            %this.corners = [158 122 85; 96 126 88; 156 116 27; 101 113 28]; % PPG
 			%this.testObj_ = Fung2013('coords', this.corners, 'plotmore', true, 'iterations', 1000, 'smoothFactor', 0.1, 'BBBuf', [10 10 4]);
 			this.testObj_ = Fung2013('coords', this.corners, 'plotmore', true, 'iterations', 100, 'BBBuf', [4 4 4]);
            this.t1w = mlfourd.ImagingContext2(fullfile(this.anatPath, 'sub-108293_20210218081030_T1w.nii.gz'));
            this.ho = mlfourd.ImagingContext2('sub-108293_20210421134537_Water_Dynamic_13_on_T1w.nii.gz');
            this.ho_sumt = mlfourd.ImagingContext2('sub-108293_20210421134537_Water_Static_12_on_T1w.nii.gz');
            %this.ho = mlfourd.ImagingContext2('hodt20190523120249_on_T1001.nii.gz');
            %this.cacheMat = fullfile(getenv('HOME'), 'MATLAB-Drive', 'mlvg', 'test', '+mlvg_unittest', 'Test_Fung2013_20210512.mat');
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
            this.testObj.buildCorners(this.corners);
            this.testObj.buildSegmentation();
            this.testObj.buildCenterlines()                
            this.testObj.buildRegistrationTargets(this.ho)
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

