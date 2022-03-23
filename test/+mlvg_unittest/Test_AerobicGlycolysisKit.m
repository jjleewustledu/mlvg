classdef Test_AerobicGlycolysisKit < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 03-Feb-2022 19:45:22 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/test/+mlvg_unittest.
    %  Developed on Matlab 9.11.0.1837725 (R2021b) Update 2 for MACI64.  Copyright 2022 John J. Lee.
    
    properties
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mlvg.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
    end
    
    methods (TestClassSetup)
        function setupAerobicGlycolysisKit(this)
            import mlvg.*
            this.testObj_ = AerobicGlycolysisKit();
        end
    end
    
    methods (TestMethodSetup)
        function setupAerobicGlycolysisKitTest(this)
            this.testObj = this.testObj_;
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
