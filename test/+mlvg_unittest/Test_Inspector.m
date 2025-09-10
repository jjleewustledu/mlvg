classdef Test_Inspector < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 08-Sep-2025 02:33:13 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/test/+mlvg_unittest.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.
    
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

        function test_plot_(this)
            %% Example usage with synthetic data

            insp = mlvg.Inspector();
            
            % Generate example data
            rng(42);  % For reproducibility
            L = 50;

            % Age vector
            a = 45 + 15 * randn(L, 1);

            % Data vectors with some correlation to age
            w = 100 + 0.5 * a + 10 * randn(L, 1);
            x = 80 - 0.3 * a + 8 * randn(L, 1);
            y = 90 + 0.2 * a + 12 * randn(L, 1);
            z = 110 + 0.7 * a + 15 * randn(L, 1);

            % Insert some NaN values randomly
            nan_indices = randperm(L, 5);
            w(nan_indices(1:2)) = NaN;
            x(nan_indices(2:3)) = NaN;
            y(nan_indices(3:4)) = NaN;
            z(nan_indices(4:5)) = NaN;

            % Create the plot
            % fig = insp.plot_age_vs_measurements_combined(a, w, x, y, z);

            % If you have your own viridis function:
            dstruct.wb = w;
            dstruct.gm = x;
            dstruct.wm = y;
            dstruct.subcortex = z;
            fig = insp.plot_age_vs_measurements_combined(a, dstruct, @viridis);
        end

        function test_reshape_to_parc(this)
            
        end
    end
    
    methods (TestClassSetup)
        function setupInspector(this)
            import mlvg.*
            this.testObj_ = Inspector();
        end
    end
    
    methods (TestMethodSetup)
        function setupInspectorTest(this)
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
