nsclassdef Test_AerobicGlycolysisKit < matlab.unittest.TestCase
    %% Tests all quantitative inferences on sub-108293/ses-20210421
    %  
    %  Created 31-Jan-2023 22:38:05 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/test/+mlvg_unittest.
    %  Developed on Matlab 9.13.0.2126072 (R2022b) Update 3 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        ho_dyn_fqfn
        immediator
        oc_avgt_fqfn
        petPath
        sessionPath
        subjectPath
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mlvg.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
        function test_quadraticKit(this)
        end
        function test_dispersiveKit(this)
        end
        function test_construct(this)
           these = mlvg.QuadraticAerobicGlycolysisKit.construct('cmrglc-posthoc');
           disp(these)
        end
        function test_constructData(this)
            theData = mlvg.QuadraticAerobicGlycolysisKit.constructData( ...
                'subjectsExpr', 'sub-108293', ...
                'tracer', 'ho', ...
                'metric', 'fs', ...
                'region', 'voxels'); 
            this.verifyEqual(theData.scanPath, ...
                '/Users/jjlee/Singularity/CCIR_01211/derivatives/sub-108293/ses-20210421/pet')
            this.verifyEqual(theData.imagingAtlas.fqfn, ...
                '/usr/local/fsl/data/standard/MNI_T1_1mm.nii.gz')
            this.verifyEqual(theData.imagingContext.fqfn, ...
                '/Users/jjlee/Singularity/CCIR_01211/derivatives/sub-108293/ses-20210421/pet/sub-108293_ses-20210421152358_trc-ho_proc-dyn_pet_on_T1w.nii.gz')
            this.verifyEqual(theData.bids.t1w_ic.fqfn, ...    
                '/Users/jjlee/Singularity/CCIR_01211/derivatives/sub-108293/ses-20210218/anat/sub-108293_ses-20210218081506_T1w_MPR_vNav_4e_RMS_orient-std.nii.gz')
            %disp(theData)
            %disp(theData.tracerOnAtlas())
        end
        function test_mlan_sessionData(this)
            setenv('SUBJECTS_DIR', '~/Singularity/CCIR_01211/derivatives')
            sd = mlan.QuadraticAerobicGlycolysisKit.constructSessionData( ...
                'cbf', ...
                'subjectsExpr', 'sub-S03292', ...
                'tracer', 'ho', ...
                'debug', true, ...
                'region', 'voxels'); 
            disp(sd)
        end
        function test_BiographData(this)
            b = mlsiemens.BiographData.createFromSession(this.immediator);
            this.verifyEqual(b.imagingContext.fileprefix, 'sub-108293_ses-20210421152358_trc-ho_proc-dyn_pet_on_T1w_decayUncorrect0')
            this.verifyEqual(seconds(b.datetimeWindow), seconds(duration('00:05:29')), RelTol=1e-4)
            this.verifyEqual(b.timeWindow, 329, RelTol=1e-4)
            %disp(b)
        end
        function test_BiographDevice(this)
            rm = mlpet.CCIRRadMeasurements.createFromSession(this.immediator);
            bd = mlsiemens.BiographVisionDevice.createFromSession( ...
                this.immediator, radMeasurements=rm);
            this.verifyEqual(bd.invEfficiency, 1.087017298092194, RelTol=1e-4)
            this.verifyEqual(bd.imagingContext.fileprefix, ...
                'sub-108293_ses-20210421152358_trc-ho_proc-dyn_pet_on_T1w_decayUncorrect0')
            plot(bd)
            %disp(bd)
        end
        function test_Ccir1211Mediator(this)
            ic = mlfourd.ImagingContext2(this.ho_dyn_fqfn);
            m = mlvg.Ccir1211Mediator(ic);

            % imagingContext
            this.verifyEqual(m.imagingContext == ic, false) % m must retain a copy of ic
            this.verifyEqual(m.imagingContext.fqfn, this.ho_dyn_fqfn)
            this.verifyEqual(m.imagingContext.size(), [208 300 320 43])

            % properties
            this.verifyEqual(m.scanPath, ic.filepath)
            this.verifyEqual(m.tracer, 'ho')

            % methods
            this.verifyEqual(m.datetime(), datetime(2021,4,21,15,23,58, TimeZone="local"));
        end
        function test_CCIRRadMeasurements(this)
            c = mlpet.CCIRRadMeasurements.createFromSession(this.immediator);
            mMR_ = table( ...
                datetime('21-Apr-2021 17:08:23', TimeZone='local'), 65.462, 28.928, 747.2, 167347, 10, 266.092, 65.462, 0.75685829041635);
            mMR_.Properties.VariableNames = c.mMR.Properties.VariableNames;
            this.verifyEqual(c.mMR{end,2:end}, mMR_{end,2:end}, RelTol=1e-4)
            this.verifyEqual(seconds(c.mMR{end,1} - mMR_{end,1}), 0, AbsTol=1e-4)
            this.verifyEqual(c.laboratory{'Hct','measurement'}, 46.8, RelTol=1e-4)
            this.verifyEqual(c.laboratory{'Hct','TimeOfMeasurement'}, ...
                datetime('21-Apr-2021 10:09:00', TimeZone='local'))
            %disp(c)
        end
        function test_immediator(this)
            this.verifyEqual(this.immediator.imagingContext.fqfn, this.ho_dyn_fqfn);
            this.verifyEqual(this.immediator.petPointSpread(), 3.566666666666666, RelTol=1e-4)
            this.verifyEqual(this.immediator.metricOnAtlas('cbf', 'voxels').fileprefix, ...
                'sub-108293_ses-20210421152358_cbf_proc-dyn_pet_on_T1w_voxels')

            % check *OnAtlas
            kit = mlpet.AbstractAerobicGlycolysisKit2.create(this.immediator);
            this.verifyEqual(kit.metricOnAtlas('cbf', tags='voxels').fileprefix, ...
                'sub-108293_ses-20210421152358_cbf_proc-dyn_pet_on_T1w_voxels')
            this.verifyEqual(kit.cbfOnAtlas(tags='voxels').fileprefix, ...
                'sub-108293_ses-20210421152358_cbf_proc-dyn_pet_on_T1w_voxels')
        end
        function test_studyRegistries(this)
            sd = mlvg.StudyData();
            disp(sd)

            reg = mlvg.Ccir1211Registry.instance();
            disp(reg)

            reg = mlan.Ccir993Registry.instance();
            disp(reg)

            reg = mlraichle.StudyRegistry.instance();
            disp(reg)
        end
        function test_TwiliteData(this)
            t = mlswisstrace.TwiliteData.createFromSession(this.immediator);
            this.verifyEqual(seconds(t.datetimeWindow), seconds(duration('00:02:24')), RelTol=1e-4)
            this.verifyEqual(t.timeWindow, 144, RelTol=1e-4)
            %disp(t)
        end
        function test_TwiliteDevice(this)
            td = mlswisstrace.TwiliteDevice.createFromSession(this.immediator);
            this.verifyClass(td.catheterKit, 'mlswisstrace.Catheter_DT20190930')
            this.verifyEqual(td.invEfficiency, 1.629152501398188, RelTol=1e-4)
            plotall(td)
            %disp(td)
        end

        function test_view(this)
            t1w = this.immediator.t1w_ic;
            wmparc = this.immediator.wmparc_on_t1w_ic;
            mask = wmparc.binarized();
            mask = mask.blurred(3.45);
            mask = mask.thresh(0.1);

            ics = containers.Map;
            ics('cbv') = mlfourd.ImagingContext2( ...
                fullfile(this.petPath, ...
                'sub-108293_ses-20210421144815_cbv_proc-dyn_pet_on_T1w_voxels.nii.gz'));
            ics('cbf') = mlfourd.ImagingContext2( ...
                fullfile(this.petPath, ...
                'sub-108293_ses-20210421152358_cbf_proc-dyn_pet_on_T1w_voxels.nii.gz'));
            ics('cmro2') = mlfourd.ImagingContext2( ...
                fullfile(this.petPath, ...
                'sub-108293_ses-20210421150523_cmro2-umol_proc-dyn_pet_on_T1w_voxels.nii.gz'));
            ics('cmrglc') = mlfourd.ImagingContext2( ...
                fullfile(this.petPath, ...
                'sub-108293_ses-20210421155638_cmrglc-umol_proc-dyn_pet_on_T1w_wmparc1.nii.gz'));
            ics('agi') = mlfourd.ImagingContext2( ...
                fullfile(this.petPath, ...
                'sub-108293_ses-20210421155638_agi_proc-dyn_pet_on_T1w.nii.gz'));
        
            for k = ics.keys
                ics(k{1}) = ics(k{1}).blurred(3.45);
                ics(k{1}) = ics(k{1}) .* mask;
                t1w.view(ics(k{1}))
            end
        
        end
    end
    
    methods (TestClassSetup)
        function setupAerobicGlycolysisKit(this)
%            setenv('SINGULARITY_HOME', '/home/usr/jjlee/mnt/CHPC_scratch/Singularity')
            import mlvg.*
            this.subjectPath = fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211/derivatives/sub-108293');
            this.sessionPath = fullfile(this.subjectPath, 'ses-20210421');
            this.petPath = fullfile(this.sessionPath, 'pet');
            this.testObj_ = []; % must call abstract factory's construct*() methods.
        end
    end
    
    methods (TestMethodSetup)
        function setupAerobicGlycolysisKitTest(this)
            this.testObj = this.testObj_;
            this.ho_dyn_fqfn = fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211/derivatives/sub-108293/ses-20210421/pet/sub-108293_ses-20210421152358_trc-ho_proc-dyn_pet_on_T1w.nii.gz');
            this.oc_avgt_fqfn = fullfile(getenv('SINGULARITY_HOMEd'), 'CCIR_01211/derivatives/sub-108293/ses-20210421/pet/sub-108293_ses-20210421144815_trc-oc_proc-dyn_pet_avgt_b25_on_T1w.nii.gz');
            this.immediator = mlvg.Ccir1211Mediator(this.ho_dyn_fqfn);
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
