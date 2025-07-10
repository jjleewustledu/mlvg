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
            local_cp = fullfile(this.local_deriv_pet, "centerline_on_pet.nii.gz'");
            remote_cp = fullfile(this.remote_deriv_pet, "centerline_on_pet.nii.gz'");

            % rsync derivatives/sub-*/ses-*/pet/centerline_on_pet.nii.gz to cluster
            assert(isfile(local_cp))
            cmd = sprintf( ...
                "rsync -ra %s %s:%s", local_cp, remote, remote_cp);
            system(cmd);
        end

    end

    %% HELPERS

    methods (Static)
        function flip_pet(folder, opts)
            arguments
                folder {mustBeFolder}
                opts.do_view logical = true
                opts.dry_run logical = false
                opts.do_static logical = true
                opts.do_dynamic logical = true
            end

            cd(folder)
            globbed = mglob("*sub-*_ses-*_trc-*createNifti*.nii.gz");
            globbed1 = [];
            if opts.do_static
                globbed1 = [globbed1, globbed(contains(globbed, "Static"))];
            end
            if opts.do_dynamic
                globbed1 = [globbed1, globbed(~contains(globbed, "Static"))];
            end
            for g = globbed1

                fprintf("%s: flip(%s,3)\n", stackstr(), g);
                if opts.dry_run
                    continue
                end

                % view Static first
                ic = mlfourd.ImagingContext2(g);
                fp = ic.fileprefix;
                ic = flip(ic, 3);  % large dynamic images keep only handles
                ic.fileprefix = fp;
                if opts.do_view && contains(ic.fileprefix, "Static")
                    ic.view();  % manually confirm
                end
                ic.save();  % overwrite silently
            end
        end

        function inspect_centerlines(folder, opts)
            arguments
                folder {mustBeFolder}
                opts.dry_run logical = false
            end

            cd(folder);
            globbed = mglob("**/sub-*_ses-*_trc-*_mipt.nii.gz");
            for mipt = globbed
                try
                    fprintf("%s: inspecting %s\n", stackstr(), mipt);
                    cl = fullfile(fileparts(mipt), "centerline_on_pet.nii.gz");
                    assert(isfile(cl))
                    if opts.dry_run
                        fprintf("fsleyes %s %s\n", mipt, cl);
                        continue
                    end
                    system(sprintf("fsleyes %s %s", mipt, cl));
                catch ME
                    handwarning(ME)
                end
            end
        end

        function populate_rawdata(opts)
            %% N.B. alternative dates of FDG scanning:
            %
            % >> populate_rawdata()
            % populate_rawdata: /data/nil-bluearc/vlassenko/RAW_IMAGES/PET/108014/**/108014*20220718*FDG*.nii.gz returned empty
            %                   ignore globbing date
            % populate_rawdata: /data/nil-bluearc/vlassenko/RAW_IMAGES/PET/108021/**/108021*20230323*FDG*.nii.gz returned empty
            % populate_rawdata: /data/nil-bluearc/vlassenko/RAW_IMAGES/PET/108121/**/108121*20231103*FDG*.nii.gz returned empty
            % populate_rawdata: /data/nil-bluearc/vlassenko/RAW_IMAGES/PET/108140/**/108140*20230424*FDG*.nii.gz returned empty
            % populate_rawdata: /data/nil-bluearc/vlassenko/RAW_IMAGES/PET/108153/**/108153*20231130*FDG*.nii.gz returned empty
            % populate_rawdata: /data/nil-bluearc/vlassenko/RAW_IMAGES/PET/108188/**/108188*20221221*FDG*.nii.gz returned empty
            %                   missing listmode from 20230828
            % populate_rawdata: /data/nil-bluearc/vlassenko/RAW_IMAGES/PET/108206/**/108206*20230327*FDG*.nii.gz returned empty
            %                   missing listmode from 20230914
            % populate_rawdata: /data/nil-bluearc/vlassenko/RAW_IMAGES/PET/108207/**/108207*20230403*FDG*.nii.gz returned empty
            %                   ignore globbing date
            % populate_rawdata: /data/nil-bluearc/vlassenko/RAW_IMAGES/PET/108242/**/108242*20230608*FDG*.nii.gz returned empty
            %                   missing listmode from 20221219
            % populate_rawdata: /data/nil-bluearc/vlassenko/RAW_IMAGES/PET/108243/**/108243*20221201*FDG*.nii.gz returned empty
            %                   ignore globbing date
            % populate_rawdata: /data/nil-bluearc/vlassenko/RAW_IMAGES/PET/108259/**/108259*20231016*FDG*.nii.gz returned empty
            % populate_rawdata: /data/nil-bluearc/vlassenko/RAW_IMAGES/PET/108266/**/108266*20230511*FDG*.nii.gz returned empty
            % populate_rawdata: /data/nil-bluearc/vlassenko/RAW_IMAGES/PET/108300/**/108300*20210517*FDG*.nii.gz returned empty
            % populate_rawdata: /data/nil-bluearc/vlassenko/RAW_IMAGES/PET/108322/**/108322*20240701*FDG*.nii.gz returned empty

            arguments
                opts.tracer {mustBeTextScalar} = "FDG"
                opts.dest_home {mustBeFolder} = pwd  % $SINGULARITY_HOME/CCIR_01211/rawdata
                opts.src_home {mustBeFolder} = "/data/nil-bluearc/vlassenko/RAW_IMAGES/PET"
                opts.ignore_globbing_date logical = false
            end

            cd(opts.dest_home);
            subs = mglob("sub-*");
            for s = subs

                pwd0 = pushd(s);

                % sub number
                sub = strrep(s, filesep, "");
                sub_num = extractAfter(sub, 4);

                % parse ses date
                sess = mglob("ses-*");
                sess = strrep(sess, filesep, "");
                ses = sess(1);
                ses_date = extractAfter(ses, 4);

                % destinaton
                dest = fullfile(opts.dest_home, sub, ses, "pet");
                ensuredir(dest);

                % glob Nick's BIDS NIfTI by sub, ses, tracer
                if opts.ignore_globbing_date
                    pattern_file = sprintf("%s*%s*.nii.gz", sub_num, opts.tracer);
                else
                    pattern_file = sprintf("%s*%s*%s*.nii.gz", sub_num, ses_date, opts.tracer);
                end
                pattern_glob = fullfile( ...
                    opts.src_home, ...
                    sub_num, ...
                    "**", ...
                    pattern_file);
                nii = mglob(pattern_glob);
                if isempty(nii)
                    fprintf("%s: %s returned empty\n", stackstr(), pattern_glob)
                    popd(pwd0);
                    continue
                end

                % copy nii and json to $SINGULARITY_HOME/CCIR_01211/rawdata
                for n = nii
                    copyfile(n, dest);
                    copyfile(strrep(n, ".nii.gz", ".json"), dest);
                end

                popd(pwd0);

            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
