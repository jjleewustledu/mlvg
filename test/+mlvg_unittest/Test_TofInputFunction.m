classdef Test_TofInputFunction < matlab.unittest.TestCase
	%% TEST_TOFINPUTFUNCTION 

	%  Usage:  >> results = run(mlraichle_unittest.Test_TofInputFunction)
 	%          >> result  = run(mlraichle_unittest.Test_TofInputFunction, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 22-Nov-2021 20:34:42 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlvg/test/+mlraichle_unittest.
 	%% It was developed on Matlab 9.11.0.1809720 (R2021b) Update 1 for MACI64.  Copyright 2021 John Joowon Lee.
 	
	properties
        bbBuffer
        bids
        cacheMat % assign after segmentations, centerlines, registration targets are ready
        corners
        do_view = true
        ho
        ho_sumt
        mmppix
 		registry
 		testObj
    end

    properties (Dependent)
        anatPath
        petPath
        sourceAnatPath
        sourcePetPath
        t1w
    end

    methods %% get/set
        function g = get.anatPath(this)
            g = this.bids.anatPath;
        end
        function g = get.petPath(this)
            g = this.bids.petPath;
        end
        function g = get.sourceAnatPath(this)
            g = this.bids.sourceAnatPath;
        end
        function g = get.sourcePetPath(this)
            g = this.bids.sourcePetPath;
        end
        function g = get.t1w(this)
            g = this.bids.t1w_ic;
        end
    end

	methods (Test)
		function test_afun(this)
 			import mlvg.*;
 			this.assumeEqual(1,1);
 			this.verifyEqual(1,1);
 			this.assertEqual(1,1);
        end
        function test_ctor(this)
            disp(this.testObj)
        end
        function test_anatomy(this)
            disp(this.testObj.anatomy)
            disp(this.testObj.anatomy_mask)
            if this.do_view
                this.testObj.anatomy.view(this.testObj.anatomy_mask)
            end
        end
        function test_segmentation(this)
        end
        function test_centerline(this)
        end
        function test_petGlobbed(this)
            this.testObj.petGlobbed();
            disp(this.testObj)
        end
        function test_taus(this)
            disp(this.testObj.taus)
            %taus = this.testObj.taus;
            disp(this.testObj.timesMid)
            %timesMid = this.testObj.timesMid;
            %save(fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211', 'taus.mat'), 'taus')
            %save(fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211', 'timesMid.mat'), 'timesMid')
        end
        function test_call(this)
            setenv('DEBUG', '')
            this.testObj_.segmentation_only = false;
            tbl = this.testObj_.call('tracerPatt', '*trc-*_proc-dyn*', ...
                'innerRadii', [0 0 0 0 ], ...
                'outerRadii', [2 4 8 16]);
            disp(tbl)
            setenv('DEBUG', '')

            % remove centerlines cache file after errors
        end
        function test_buildPetOnTof(this)
            import mlfourd.*
            petfile = fullfile(this.petPath, 'sub-108293_ses-20210421155709_trc-fdg_proc-dyn_pet_on_T1w.4dfp.hdr');
            petic = this.testObj_.buildPetOnTof(petfile);
            petic.view(this.testObj.anatomy);
        end
        function test_buildMasks(this)
            this.testObj.T1_mask.view(this.bids.T1_ic);
            this.testObj.t1w_mask.view(this.bids.t1w_ic);
            this.testObj.tof_mask.view(this.bids.tof_ic);
        end
        function test_flirtT1wOnTof(this)
            f = this.testObj.flirtT1wOnTof();
            disp(f)
        end
        function test_ensureBoxInFieldOfView(this)
            ic = this.testObj_.anatomy;
            pc = ic.pointCloud('threshp', 25);
            figure; pcshow(pc)
            X = round(pc.Location(:,1));
            Y = round(pc.Location(:,2));
            Z = round(pc.Location(:,3));
            Z(Z < 64) = -1;
            [X,Y,Z] = this.testObj_.ensureSubInFieldOfView(X, Y, Z);
            
            ind = sub2ind(size(ic), X, Y, Z);
            img = zeros(size(ic));
            img(ind) = 1;
            ifc = ic.nifti;
            ifc.img = img;
            ifc.view();
        end
	end

 	methods (TestClassSetup)
		function setupTofInputFunction(this)
 			import mlvg.*;
            this.bids = Ccir1211Bids.create( ...
                'destinationPath', fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211', 'derivatives', 'sub-108293', 'ses-20210421', 'pet', ''));
            this.mmppix = [0.260417 0.260417 0.5];

            cd(this.petPath) 
            this.ho = mlfourd.ImagingContext2( ...
                fullfile(this.petPath, 'sub-108293_ses-20210421152358_trc-ho_proc-dyn_pet_on_T1w.nii.gz'));
            this.ho_sumt = mlfourd.ImagingContext2( ...
                fullfile(this.petPath, 'sub-108293_ses-20210421152358_trc-ho_proc-static_pet_on_T1w.nii.gz'));

            % PPG T1001; [ x y z; ... ]; [ [RS medial]; [LS medial]; [RI lateral]; [LI lateral] ].  NIfTI space.
            this.corners = [354 511 87; 309 511 87; 467 398 36; 201 398 36];
            this.bbBuffer = ceil([0 0 0] ./ this.mmppix);

            this.project_corners_ = containers.Map; % voxel coords
            this.project_corners_('sub-108293') = this.corners;
            this.project_bbBuffer_ = containers.Map; % lengths in mm
            this.project_bbBuffer_('sub-108293') = this.bbBuffer;
            this.project_iterations_ = containers.Map;
            this.project_iterations_('sub-108293') = 50;
            this.project_seg_thresh_ = containers.Map;
            this.project_seg_thresh_('sub-108293') = 190;
            this.project_contract_bias_ = containers.Map;
            this.project_contract_bias_('sub-108293') = 0.2;

 			this.testObj_ = TofInputFunction('corners', this.corners, 'iterations', 100, 'bbBuffer', this.bbBuffer, ...
                'contractBias', 0.2, 'smoothFactor', 0, 'segmentationThresh', 190, 'segmentationOnly', false, ...
                'innerRadius', 0, 'outerRadius', 2, 'subjectFolder', 'sub-108293', 'plotdebug', true, 'plotclose', true, ...
                'destinationPath', this.petPath);
 		end
	end

 	methods (TestMethodSetup)
		function setupTofInputFunctionTest(this)
            if isempty(this.testObj_)
                return
            end
 			this.addTeardown(@this.cleanTestMethod);
            this.testObj = copy(this.testObj_);
 		end
	end

	properties (Access = private)
        project_bbBuffer_
        project_corners_
        project_iterations_
        project_seg_thresh_
        project_contract_bias_
 		testObj_
 	end

	methods (Access = private)
		function cleanTestMethod(this)
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

