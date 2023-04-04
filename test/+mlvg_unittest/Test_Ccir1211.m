classdef Test_Ccir1211 < matlab.unittest.TestCase
    %% See also:  mlvg_unittest.Test_PETDirector, mlvg_unittest.Test_PETBuilder.
    %  
    %  Created 13-Jun-2022 18:56:48 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/test/+mlvg_unittest.
    %  Developed on Matlab 9.12.0.1956245 (R2022a) Update 2 for MACI64.  Copyright 2022 John J. Lee.
    
    properties
        sesObj
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mlvg.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
        function test_setup(this)
            disp(this)
            disp(this.sesObj)
            disp(this.testObj.bids)
            disp(this.testObj.study)
        end
        function test_minimal_session(this)
            sesObj_ = mlvg.SessionData( ...
                'subjectData', mlvg.SubjectData('subjectFolder', 'sub-108293'), ...
                'sessionFolder', 'ses-20210421', ...
                'scanIndex', 1, ...
                'tracer', 'HO', ...
                'ac', true, ...
                'region', 'voxels', ...
                'metric', 'fs');
            disp(sesObj_)
        end

        %% mlvg.Ccir1211

        function test_call(this)
        end
        function test_call_sub(this)
        end
        function test_call_ses(this)
        end
        function test_registry(this)
            this.verifyNotEmpty(mlvg.Ccir1211Registry.instance())
        end

        %% mlkinetics.OxygenMetabKit

        function test_OxygenMetabKit(this)
        end        
        function test_oxy_make_tracer_iterator(this)
        end
        function test_oxy_make_scanner(this)
        end
        function test_oxy_make_input_function(this)
        end
        function test_oxy_make_model(this)
        end

        %% mlkinetics.GlucoseMetabKit

        function test_GlucoseMetabKit(this)
        end        
        function test_glc_make_tracer(this)
        end
        function test_glc_make_scanner(this)
        end
        function test_glc_make_input_function(this)
        end
        function test_glc_make_model(this)
        end
    end
    
    methods (TestClassSetup)
        function setupCcir1211(this)
            import mlvg.*
            this.testObj_ = Ccir1211();
        end
    end
    
    methods (TestMethodSetup)
        function setupCcir1211Test(this)
            this.testObj = this.testObj_;
            this.sesObj = mlvg.SessionData( ...
                'subjectData', mlvg.SubjectData('subjectFolder', 'sub-108293'), ...
                'sessionFolder', 'ses-20210421', ...
                'scanIndex', 1, ...
                'tracer', 'HO', ...
                'ac', true, ...
                'region', 'voxels', ...
                'metric', 'fs');  
            this.addTeardown(@this.cleanTestMethod)
        end
    end
    
    properties (Access = private)
        testObj_
    end
    
    methods (Access = private)
        function cleanTestMethod(this)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
