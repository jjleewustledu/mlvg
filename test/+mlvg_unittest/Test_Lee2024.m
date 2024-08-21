classdef Test_Lee2024 < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 23-Jul-2024 23:07:07 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/test/+mlvg_unittest.
    %  Developed on Matlab 24.1.0.2653294 (R2024a) Update 5 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        testObj
    end

    methods (Static)
        function s = mov2static(s)
            s = strrep(s, "createNiftiMovingAvgFrames", "createNiftiStatic");
        end
    end
    
    methods (Test)
        function test_afun(this)
            import mlvg.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
        function test_build_maframes_timeAppend_parc(this)
            %% build schaeffer representations for all the CO dynamic data, all frame window lengths
            %  assumes filesystem has canonical pet folder with tau-tags in filenames, compatible with BidsKit &
            %  ImagingMediator.schaeffer_ic

            lee2024 = this.testObj;

            taus = [10,5,3];
            parfor tidx = 1:length(taus)
                globbed = asrow(glob(sprintf( ...
                    '/home/usr/jjlee/mnt/Singularity/CCIR_01211/sourcedata/sub-*/ses-*/pet/sub-*_ses-*_trc-co*-tau%i-createNiftiMovingAvgFrames.nii.gz', taus(tidx))));
                for gidx = 1:length(globbed)
                    try
                        % pwd0 = pushd(fileparts(fileparts(fileparts(g{1}))));  % sub-*/ses-*/                        
                        lee2024.build_maframes_timeAppend_parc(globbed{gidx}); %#ok<PFBNS>
                        % popd(pwd0);
                    catch ME
                        handwarning(ME)
                    end
                end
            end
        end
        function test_build_maframes_timeAppend_parc_sub108237(this)
            %% build schaeffer representations for sub-108237 only; unit test

            pth = "/home/usr/jjlee/mnt/Singularity/CCIR_01211/sourcedata/sub-108237/ses-20221031100910/pet";
            fn = "sub-108237_ses-20221031100910_trc-co_proc-delay0-BrainMoCo2-tau10-createNiftiMovingAvgFrames.nii.gz";
            fqfn = fullfile(pth, fn);

            ic = this.testObj.build_maframes_timeAppend_parc(fqfn);
            ic.view();
        end
        function test_build_maframes_timeAppend_parc_ori(this)
            %% build schaeffer representations for all the CO dynamic data, all frame window lengths
            %  reorganize filesystem that has folder tau=*s, which are incompatible with BidsKit &
            %  ImagingMediator.schaeffer_ic

            for tau = [10,5,3]
                globbed = asrow(glob(sprintf( ...
                    '/home/usr/jjlee/mnt/Singularity/CCIR_01211/sourcedata/sub-*/ses-*/tau=%is/pet/sub-*_ses-*_trc-co*createNiftiMovingAvgFrames.nii.gz', tau)));
                for g = globbed
                    try
                        pwd0 = pushd(fileparts(fileparts(fileparts(g{1}))));  % sub-*/ses-*/
                        ensuredir("pet");
                        cd("pet");
                        fqfp = myfileprefix(g{1});
                        fp = mybasename(g{1});
                        fp1 = strrep(fp, "-create", "-tau"+tau+"-create");
                        if isfile(fqfp+".nii.gz")
                            system(sprintf("mv -f %s.nii.gz %s.nii.gz", fqfp, fp1));
                            system(sprintf("mv -f %s.json %s.json", fqfp, fp1));
                        end
                        if isfile(this.mov2static(fqfp)+".nii.gz")
                            system(sprintf("mv -f %s.nii.gz %s.nii.gz", this.mov2static(fqfp), this.mov2static(fp1)));
                            system(sprintf("mv -f %s.json %s.json", this.mov2static(fqfp), this.mov2static(fp1)));
                        end
                        this.testObj.build_maframes_timeAppend_parc(fullfile(pwd, fp1+".nii.gz"));
                        popd(pwd0);
                    catch ME
                        handwarning(ME)
                    end
                end
            end
        end
    end
    
    methods (TestClassSetup)
        function setupLee2024(this)
            import mlvg.*
            this.testObj_ = Lee2024();
        end
    end
    
    methods (TestMethodSetup)
        function setupLee2024Test(this)
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
