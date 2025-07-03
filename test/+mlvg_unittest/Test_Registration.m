classdef Test_Registration < matlab.unittest.TestCase
	%% TEST_REGISTRATION 

	%  Usage:  >> results = run(mlvg_unittest.Test_Registration)
 	%          >> result  = run(mlvg_unittest.Test_Registration, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 10-Nov-2021 14:12:37 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlvg/test/+mlvg_unittest.
 	%% It was developed on Matlab 9.11.0.1769968 (R2021b) for MACI64.  Copyright 2021 John Joowon Lee.
 	
	properties
        bids
 		testObj
 	end

	methods (Test)
		function test_afun(this)
 			import mlvg.*;
 			this.assumeEqual(1,1);
 			this.verifyEqual(1,1);
 			this.assertEqual(1,1);
        end
        function test_ctor(this)
            disp(this.bids)
            disp(this.testObj)
            o = this.testObj;
            o.t1w_ic.view(o.t1w_weight_ic, o.wmparc_on_t1w_ic)
        end
        function test_flirt_t1w_to_tof(this)
            o = this.testObj;
            ic = o.flirt_t1w_to_tof();
            ic.view(this.bids.tof_ic.fqfn)
        end
        function test_flirt_tof_to_t1w(this)
            o = this.testObj;
            ic = o.flirt_tof_to_t1w();
            ic.view(this.bids.t1w_ic.fqfn)
        end
        function test_flirt_all_to_t1w(this)
            this = mlvg.Registration.flirt_all_to_t1w( ...
                "projectPath", fullfile(getenv("SINGULARITY_HOME"), "CCIR_01211"), ...
                "subjectFolder", "sub-108007");
            pwd0 = pushd(this.bids.anatPath);
            this.t1w_ic.view( ...
                "sub-108007_20201203150706_T2w_on_T1w.nii.gz", ...
                "sub-108007_20201203145958_FLAIR_on_T1w.nii.gz", ...
                "sub-108007_20201203154051_angio_on_T1w.nii.gz")
            popd(pwd0)
        end
	end

 	methods (TestClassSetup)
		function setupRegistration(this)
 			import mlvg.*;
            cd(fullfile(getenv("SINGULARITY_HOME"), "CCIR_01211", "derivatives", "sub-108007", "ses-20201203", "anat", ""));
            this.bids = mlvg.Ccir1211Bids( ...
                "projectPath", fullfile(getenv("SINGULARITY_HOME"), "CCIR_01211"), ...
                "subjectFolder", "sub-108007");
 			this.testObj_ = Registration(this.bids);
 		end
	end

 	methods (TestMethodSetup)
		function setupRegistrationTest(this)
 			this.testObj = this.testObj_;
 			this.addTeardown(@this.cleanTestMethod);
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

