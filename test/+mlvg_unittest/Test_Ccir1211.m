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
        function test_ccir1211mediator(this)
            fqfn = fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211/derivatives/sub-108293/ses-20210421/pet', ...
                'sub-108293_ses-20210421171325_trc-fdg_proc-static-phantom_pet.nii.gz');
            med = mlvg.Ccir1211Mediator(fqfn);
            
            this.verifyClass(med, 'mlvg.Ccir1211Mediator')
            this.verifyEqual(med.t1w_ic.fileprefix, 'sub-108293_ses-20210218081506_T1w_MPR_vNav_4e_RMS_orient-std')
            this.verifyEqual(med.tof_ic.fileprefix, 'sub-108293_ses-20210218092914_tof_fl3d_tra_p2_multi-slab_orient-std')
            this.verifyEqual(med.imagingAtlas.fileprefix, 'MNI152_T1_1mm')            

            this.verifyEqual(med.imagingContext.filepath, ...
                fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211/derivatives/sub-108293/ses-20210421/pet'))
            this.verifyEqual(med.imagingContext.fileprefix, ...
                'sub-108293_ses-20210421171325_trc-fdg_proc-static-phantom_pet')
            this.verifyEqual(med.timeOffsetConsole, seconds(31))
        end
        function test_derivatives(this)
            fqfn = fullfile(getenv("SINGULARITY_HOME"), ...
                "CCIR_01211", "sourcedata", "sub-108293", 'ses-20210421152358', "pet", ...
                "sub-108293_ses-20210421152358_trc-ho_proc-for-testing.nii.gz");
            med = mlvg.Ccir1211Mediator(fqfn);
            icd = med.prepare_derivatives();
            this.verifyTrue(isfile(icd.fqfn));
            deleteExisting(icd.fqfn);
        end
        function test_datetime(this)
            fqfn = fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211/derivatives/sub-108293/ses-20210421/pet', ...
                'sub-108293_ses-20210421171325_trc-fdg_proc-static-phantom_pet.nii.gz');
            med = mlvg.Ccir1211Mediator(fqfn);
            
            this.verifyEqual(med.timeOffsetConsole, seconds(31))
            this.verifyEqual(datetime(med), datetime(2021,04,21, 17,13,25, TimeZone='local'))
            this.verifyEqual(datetime_console_adjusted(med), datetime(2021,04,21, 17,12,54, TimeZone='local'))
        end
        function test_initialize(this)
            fqfn = fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211/derivatives/sub-108293/ses-20210421/pet', ...
                'sub-108293_ses-20210421171325_trc-fdg_proc-static-phantom_pet.nii.gz');
            med = mlvg.Ccir1211Mediator(fqfn);

            fqfn1 = fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211/derivatives/sub-108293/ses-20210421/pet', ...
                'sub-108293_ses-20210421155709_trc-fdg_proc-dyn_pet_keepframes-4-5_avgt_b25.nii.gz');
            med.initialize(fqfn1);
            %disp(med)
            %disp(med.imagingContext)
            this.verifyEqual(med.imagingContext.fileprefix, ...
                'sub-108293_ses-20210421155709_trc-fdg_proc-dyn_pet_keepframes-4-5_avgt_b25')
        end
        function test_minimal_legacy_session(this)
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
            fdg_fqfn = fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211/derivatives/sub-108293/ses-20210421/pet', ...
                'sub-108293_ses-20210421155709_trc-fdg_proc-dyn_pet.nii.gz');
            kit = mlkinetics.GlcMetabKit.create_cmrglc( ...
                bids_tags='ccir1211', bids_fqfn=fdg_fqfn);
        end  
        function test_glc_make_bids(this)
        end
        function test_glc_make_input_function(this)
        end
        function test_glc_make_model(this)
        end
        function test_glc_make_parc(this)
        end
        function test_glc_make_scanner(this)
        end
        function test_glc_make_tracer(this)
        end
    end
    
    methods (TestClassSetup)
        function setupCcir1211(this)
        end
    end
    
    methods (TestMethodSetup)
        function setupCcir1211Test(this)
            % this.testObj = mlvg.Ccir1211;
            % this.sesObj = mlvg.SessionData( ...
            %     'subjectData', mlvg.SubjectData('subjectFolder', 'sub-108293'), ...
            %     'sessionFolder', 'ses-20210421', ...
            %     'scanIndex', 1, ...
            %     'tracer', 'HO', ...
            %     'ac', true, ...
            %     'region', 'voxels', ...
            %     'metric', 'fs');  
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
