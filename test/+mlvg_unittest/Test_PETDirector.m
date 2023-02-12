classdef Test_PETDirector < matlab.unittest.TestCase
    %% See also:  mlvg_unittest.Test_PETBuilder, mlvg_unittest.Test_Ccir1211.
    %  
    %  Created 05-Sep-2022 15:50:41 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/test/+mlvg_unittest.
    %  Developed on Matlab 9.12.0.2039608 (R2022a) Update 5 for MACI64.  Copyright 2022 John J. Lee.
    
    properties
        study_home
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mlvg.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
        function test_announce(~)
            announce('begin')
            t = tic;
            announce('end', t)
        end
        function test_ctor(this)
            pd = mlvg.PETDirector('workdir', fullfile(this.study_home, 'rawdata', 'sub-108287', ''));
            self_test(pd);
            pd = mlvg.PETDirector('workdir', fullfile(this.study_home, 'sourcedata', 'sub-108287', ''));
            self_test(pd);
            pd = mlvg.PETDirector('workdir', fullfile(this.study_home, 'derivatives', 'sub-108287', ''));
            self_test(pd);
        end

        function test_build_study(this)
            this.testObj.build_study();
        end
        function test_build_filesystem(this)
            this.testObj.build_filesystem();
        end
        function test_build_unpacked(this)
            this.testObj.build_unpacked();
        end
        function test_build_reconstruction_nipet(this)
            obj = mlvg.PetDirector('workdir',  ...
                '/data/anlab/jjlee/Singularity/CCIR_00993/derivatives/nipet/ses-E268533');
            obj.build_reconstruction();
        end
        function test_build_reconstruction_jsrecon(this)
            obj = mlvg.PETDirector('workdir',  ...
                fullfile(this.study_home, 'rawdata', 'sub-108007'));
            disp(obj)
            disp(obj.reconstruction_builder);
            obj.build_reconstruction();
        end

        function test_build_bids(this)
            this.testObj.build_bids();
        end
        function test_build_session(this)
            b = this.testObj.build_session();
            this.verifyClass(b, 'mlpipeline.ISessionBuilder')
        end
        function test_build_subject(this)
            b = this.testObj.build_subject();
            this.verifyClass(b, 'mlpipeline.ISubjectBuilder')
        end
        
        function test_build_calibration(this)
            b = this.testObj.build_calibration();
            this.verifyClass(b, 'mlpipeline.ICalibrationBuilder')
        end
        function test_build_tracer(this)
            this.testObj.build_tracer();
        end
        function test_build_imaging(this)
            b = this.testObj.build_imaging();
            this.verifyClass(b, 'mlpipeline.IImagingBuilder')
        end
        function test_build_input_function(this)
            b = this.testObj.build_input_functiong();
            this.verifyClass(b, 'mlpipeline.IInputFunctionBuilder')
        end

        function test_build_alignment(this)
            b = this.testObj.build_alignment();
            this.verifyClass(b, 'mlpipeline.IAlignmentBuilder')
        end
        function test_build_model(this)
            b = this.testObj.build_model();
            this.verifyClass(b, 'mlpipeline.IModelBuilder')
        end
        function test_build_report(this)
            b = this.testObj.build_report();
            this.verifyClass(b, 'mlpipeline.IReportBuilder')
        end
        function test_build_quality_assurance(this)
            b = this.testObj.build_quality_assurance();
            this.verifyClass(b, 'mlpipeline.IQualityAssuranceBuilder')
        end
    end
    
    methods (TestClassSetup)
        function setupPETDirector(this)
            import mlvg.*
            this.study_home = fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211');
            this.testObj_ = PETDirector( ...
                'workdir', fullfile(this.study_home, 'rawdata', 'sub-108187', 'cnda.wustl.edu/108187_PET2_20220822_FDG'), ...
                'reconstruction_builder', 'jsrecon');
        end
    end
    
    methods (TestMethodSetup)
        function setupPETDirectorTest(this)
            this.testObj = this.testObj_;
            self_test(this.testObj);
            %this.addTeardown(@this.cleanTestMethod)
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
