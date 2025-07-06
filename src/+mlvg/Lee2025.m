classdef Lee2025 < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 04-Jul-2025 01:08:28 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.
    

    properties
        nii
        out_dir
    end

    properties (Dependent)
        local_deriv_pet
        local_nmaf
        local_out_dir
        local_src_pet
        remote_deriv_pet
    end

    methods  %% GET
        function g = get.local_deriv_pet(this)
            nii_src_pet = myfileparts(this.nii);
            nii_deriv_pth = strrep(nii_src_pet, "sourcedata", "derivatives");
            g = fullfile(this.local_out_dir, nii_deriv_pth);
            ensuredir(g);
        end

        function g = get.local_nmaf(this)
            [nii_src_pth,nii_fp] = myfileparts(this.nii);
            g = fullfile(this.local_out_dir, nii_src_pth, nii_fp + ".nii.gz");
        end

        function g = get.local_out_dir(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "CCIR_01211");
            ensuredir(g);
        end

        function g = get.local_src_pet(this)
            nii_src_pet = myfileparts(this.nii);
            g = fullfile(this.local_out_dir, nii_src_pet);
            ensuredir(g)
        end

        function g = get.remote_deriv_pet(this)
            nii_src_pet = myfileparts(this.nii);
            nii_deriv_pet = strrep(nii_src_pet, "sourcedata", "derivatives");
            g = fullfile(this.out_dir, nii_deriv_pet);
        end
    end

    methods
        function this = Lee2025(nii, opts)
            arguments
                nii {mustBeTextScalar} = "sourcedata/sub-108007/ses-20210219145054/pet/sub-108007_ses-20210219145054_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz"
                opts.out_dir {mustBeText} = "/scratch/jjlee/Singularity/CCIR_01211"
            end
            this.nii = nii;
            this.out_dir = opts.out_dir;
            % assert(isfile(fullfile(this.out_dir, this.nii)))
        end

        function this = call_ifk(this, opts)
            arguments
                this mlvg.Lee2025
                opts.nii = this.nii;
                opts.out_dir = this.out_dir;
                opts.method {mustBeTextScalar} = "do_make_input_func"
                opts.reference_tracer {mustBeTextScalar} = "ho"
                opts.steps {mustBeNumericOrLogical} = 1
            end
            lsteps = false(1, 6);
            lsteps(opts.steps) = true;

            % generate bids_fqfn
            bids_fqfn = fullfile(opts.out_dir, opts.nii);

            % generate ifk
            bk = mlkinetics.BidsKit.create( ...
                bids_tags="ccir1211", ...
                bids_fqfn=bids_fqfn);
            tk = mlkinetics.TracerKit.create( ...
                bids_kit=bk, ...
                ref_source_props=datetime(2022,2,1, TimeZone="local"), ...
                tracer_tags="", ...
                counter_tags="caprac");
            sk = mlkinetics.ScannerKit.create( ...
                bids_kit=bk, ...
                tracer_kit=tk, ...
                scanner_tags="vision");
            ifk = mlkinetics.InputFuncKit.create( ...
                bids_kit=bk, ...
                tracer_kit=tk, ...
                scanner_kit=sk, ...
                input_func_tags="mipidif", ...
                input_func_fqfn="");

            delete_large_files = any(opts.steps == length(lsteps));
            ifk.(opts.method)( ...
                steps=lsteps, delete_large_files=delete_large_files, reference_tracer=opts.reference_tracer);
        end

        function this = draw_centerline(this)
            this.pull_for_centerline();
            this = this.call_ifk( ...
                nii=this.local_nmaf, ...
                out_dir=this.local_out_dir, ...
                steps=2);
            this.push_new_centerline();
        end

        function this = pull_for_centerline(this)
            [nii_src_pet,nii_fp] = myfileparts(this.nii);
            remote = "login3.chpc.wustl.edu";
            remote_nmaf = fullfile(this.out_dir, nii_src_pet, nii_fp + ".*");            

            % rsync *createNiftiMovingAvgFrames.nii.gz to local
            if ~isfile(this.local_nmaf)
                cmd = sprintf( ...
                    "rsync -a %s:%s %s", remote, remote_nmaf, this.local_src_pet);
                system(cmd);
                assert(isfile(this.local_nmaf))
            end
            
            % rsync derivatives/sub-*/ses-*/pet to local
            if ~isfile(fullfile(this.local_deriv_pet, nii_fp + "_mipt.nii.gz"))
                cmd = sprintf( ...
                    "rsync -ra %s:%s/ %s/", remote, this.remote_deriv_pet, this.local_deriv_pet);
                system(cmd);
                assert(isfile(fullfile(this.local_deriv_pet, nii_fp + "_mipt.nii.gz")))
            end
        end 

        function this = push_new_centerline(this)
            local_cp = fullfile(this.local_deriv_pet, "centerline_pet.nii.gz'");
            remote_cp = fullfile(this.remote_deriv_pet, "centerline_pet.nii.gz'");

            % rsync derivatives/sub-*/ses-*/pet/centerline_pet.nii.gz to cluster
            assert(isfile(local_cp))
            cmd = sprintf( ...
                "rsync -ra %s %s:%s", local_cp, remote, remote_cp);
            system(cmd);
        end

    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
