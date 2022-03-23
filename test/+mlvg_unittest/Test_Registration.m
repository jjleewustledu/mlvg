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
 		registry
 		testObj
 	end

	methods (Test)
		function test_afun(this)
 			import mlvg.*;
 			this.assumeEqual(1,1);
 			this.verifyEqual(1,1);
 			this.assertEqual(1,1);
 		end
	end

 	methods (TestClassSetup)
		function setupRegistration(this)
 			import mlvg.*;
            cd(fullfile(getenv("SINGULARITY_HOME"), "CCIR_01211", "derivatives", "sub-108293", "ses-20210218", "anat", ""));
            this.bids = mlvg.Ccir1211Bids( ...
                "projectPath", fullfile(getenv("SINGULARITY_HOME"), "CCIR_01211"), ...
                "subjectFolder", "sub-108293");
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

