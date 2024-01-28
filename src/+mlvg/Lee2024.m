classdef Lee2024 < handle & mlsystem.IHandle
    %% For manuscript to Physics in Biology & Medicine, 2024.
    %  
    %  Created 03-Jan-2024 21:45:43 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        
    end

    methods
        function draw_centerlines(this, row)            
            T = this.table_maframes_nii();
            pth = myfileparts(T.bids_fqfn(row));
            pth = strrep(pth, "sourcedata", "derivatives");
            pwd0 = pushd(pth);

            lett = ["x", "y", "z"];
            for ax = 3:-1:1
                tof = mglob(sprintf("tof_on_sub*_max%g.nii.gz", ax));
                mipt = mglob(sprintf("sub*createNiftiMovingAvgFrames_mipt_max%g.nii.gz", ax));
                assert(1 == length(tof));
                assert(1 == length(mipt));
                tof_ic = mlfourd.ImagingContext2(tof);
                mipt_ic = mlfourd.ImagingContext2(mipt);

                fprintf("Please draw arterial centerlines as overlaid images on the tof and mipt.\n" + ...
                    "Save files: centerline_%sl and centerline_%sr.\n\n", lett(ax));
                tof_ic.view(mipt_ic) 
            end

            popd(pwd0);
        end
        function dilate_centerlines(this)
            T = this.table_centerlines();
            for cidx = 1:size(T, 1)
                ic = mlfourd.ImagingContext2(T.bids_fqfn{cidx});
                for m = 2:10
                    this.dilate_centerline(ic, m);
                end
            end
        end
        function ic = dilate_centerline(~, ic, m)
            assert(m > 1)

            p = floor(mean(m, m)/2);
            se = strel('cuboid', [m,m,p]);
            fp = ic.fileprefix;
            ic = ic.imdilate(se);
            ic.fileprefix = sprintf("%s_dilate-m%i", fp, m);
            ic.save();
        end
        function cl = shift_centerline(this, row, opts)
            arguments
                this mlvg.Lee2024
                row double
                opts.s double = 2
                opts.sz double = [440,440,159]
            end

            T = this.table_maframes_nii();
            pth = myfileparts(T.bids_fqfn(row));
            pth = strrep(pth, "sourcedata", "derivatives");
            pwd0 = pushd(pth);

            copyfile("centerline_on_pet.nii.gz", "centerline_on_pet_previous.nii.gz");
            cl = mlfourd.ImagingFormatContext2("centerline_on_pet.nii.gz");
            img = cl.img;
            Dx = opts.sz(1); % 440
            Dx2 = floor(opts.sz(1)/2); % 220
            cl.img(Dx2+1:Dx-opts.s,:,:) = img(Dx2+1+opts.s:Dx,:,:);
            cl = mlfourd.ImagingContext2(cl);
            cl.save();

            popd(pwd0);            
        end
        function cl = view_centerline(this, row)            
            T = this.table_maframes_nii();
            pth = myfileparts(T.bids_fqfn(row));
            pth = strrep(pth, "sourcedata", "derivatives");
            pwd0 = pushd(pth);
            
            fprintf("%s: %g\n", stackstr(), row)

            try
                cl = mlfourd.ImagingContext2("centerline_on_pet.nii.gz");
                mipt = mglob("*_mipt.nii.gz");
                if any(contains(mipt, "oo"))
                    disp(mipt(end))
                    mipt = mlfourd.ImagingContext2(mipt(end));
                    static = mglob("sub-*-createNiftiStatic.nii.gz");
                    static = mlfourd.ImagingContext2(static(end));
                    tof = mglob("tof_on_*-createNiftiStatic.nii.gz");
                    tof = mlfourd.ImagingContext2(tof(end));
                    t1w = mglob("T1w_on_*-createNiftiStatic.nii.gz");
                    t1w = mlfourd.ImagingContext2(t1w(end));
                    t1w.view(tof, static, mipt, cl)
                else
                    disp(mipt(1))
                    mipt = mlfourd.ImagingContext2(mipt(1));
                    static = mglob("sub-*-createNiftiStatic.nii.gz");
                    static = mlfourd.ImagingContext2(static(1));
                    tof = mglob("tof_on_*-createNiftiStatic.nii.gz");
                    tof = mlfourd.ImagingContext2(tof(1));
                    t1w = mglob("T1w_on_*-createNiftiStatic.nii.gz");
                    t1w = mlfourd.ImagingContext2(t1w(1));
                    t1w.view(tof, static, mipt, cl)
                end
                
            catch ME
                handwarning(ME)
            end

            popd(pwd0);
        end
        function build_all_twilite(this)
            
        end
        function call_ifk(this, row, opts)
            arguments
                this mlvg.Lee2024
                row double
                opts.steps logical = [1,1,0,0,0,0]
                opts.delete_large_files logical = false
                opts.reference_tracer {mustBeTextScalar} = "fdg"
            end

            try
                T = this.table_maframes_nii();
                bids_fqfn = T.bids_fqfn(row);
                if any(opts.steps(3:end)) % prep to use all time frames for IDIF
                    bids_fqfn = strrep(bids_fqfn, "-delay30-", "-delay0-");
                end
                % if ~any(opts.steps(5:end))
                %     bids_fqfn = strrep(bids_fqfn, "-createNiftiMovingAvgFrames", "-createNiftiStatic");
                % end
                assert(isfile(bids_fqfn)) % sanity
                deriv_pth = strrep(fileparts(bids_fqfn), "sourcedata", "derivatives");
                ensuredir(deriv_pth);
                %deleteExisting(fullfile(deriv_pth, "*proc-MipIdif*"))

                if strcmp(T.tracer(row), "fdg")
                    tracer_tags = "18F";
                else
                    tracer_tags = "15O";
                end

                ifk = mlkinetics.InputFuncKit.create_from_tags( ...
                    bids_fqfn=bids_fqfn, ...
                    bids_tags="ccir1211", ...
                    scanner_tags="vision", ...
                    tracer_tags=tracer_tags, ...
                    input_func_tags="mipidif");
                ifk.do_make_input_func( ...
                    steps=opts.steps, reference_tracer=opts.reference_tracer, delete_large_files=opts.delete_large_files);
            catch ME
                handwarning(ME)
            end
        end
        function plot_bigraph(this, trc, ti)
            Ttwil = this.table_twilite_nii();
            Tidif = this.table_idif_nii();

            Utwil = Ttwil(Ttwil.tracer == trc, :);
            Uidif = Tidif(Tidif.tracer == trc, :);

            figure; hold on
            for r = 1:size(Uidif, 1)
                plot(Uidif.timesMid{r}, Uidif.img{r}, LineWidth=2); end; hold off
            ylabel("activity (kBq/mL)")
            fontsize(scale = 1.5)
            title(ti)
            xlim(this.xlim(trc));

            figure; hold on
            for r = 1:size(Utwil, 1)
                plot(Utwil.timesMid{r}, Utwil.img{r}, LineWidth=2); end; hold off
            set(gca, 'YDir', 'reverse')
            xlabel("time (s)")
            fontsize(scale = 1.5)
            xlim(this.xlim(trc));
        end
        function plot_bigraph_nest(this, trc, ti)
            Ttwil = this.table_twilite_nest();
            Tidif = this.table_idif_nest();

            Utwil = Ttwil(Ttwil.tracer == trc, :);
            Uidif = Tidif(Tidif.tracer == trc, :);

            figure; hold on
            for r = 1:size(Uidif, 1)
                plot(Uidif.timesMid{r}, Uidif.img{r}, LineWidth=2); end; hold off
            ylabel("activity (kBq/mL)")
            fontsize(scale = 1.5)
            title(ti)
            xlim(this.xlim(trc));

            figure; hold on
            for r = 1:size(Utwil, 1)
                plot(Utwil.timesMid{r}, Utwil.img{r}, LineWidth=2); end; hold off
            set(gca, 'YDir', 'reverse')
            xlabel("time (s)")
            fontsize(scale = 1.5)
            xlim(this.xlim(trc));
        end
        function T = table_centerlines(this)
            arguments
                this mlvg.Lee2024
            end

            if ~isempty(this.table_centerlines_)
                T = this.table_centerlines_;
                return
            end

            T = this.table_maframes_nii;
            T.bids_fqfn = fullfile( ...
                myfileparts(strrep(T.bids_fqfn, "sourcedata", "derivatives")), "centerline_on_pet");
        end
        function T = table_maframes_nii(this)
            if ~isempty(this.table_maframes_nii_)
                T = this.table_maframes_nii_;
                return
            end

            console_clock(1:5) = seconds(31);
            bids_fqfn(1) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108293", "ses-20210421144815", "pet", ...
                "sub-108293_ses-20210421144815_trc-co_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(2) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108293", "ses-20210421150523", "pet", ...
                "sub-108293_ses-20210421150523_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(3) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108293", "ses-20210421152358", "pet", ...
                "sub-108293_ses-20210421152358_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(4) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108293", "ses-20210421154248", "pet", ...
                "sub-108293_ses-20210421154248_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(5) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
               "CCIR_01211", "sourcedata", "sub-108293", "ses-20210421155709", "pet", ...
               "sub-108293_ses-20210421155709_trc-fdg_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");           

            console_clock(6:10) = seconds(17);
            bids_fqfn(6) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108237", "ses-20221031100910", "pet", ...
                "sub-108237_ses-20221031100910_trc-co_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(7) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108237", "ses-20221031102320", "pet", ...
                "sub-108237_ses-20221031102320_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(8) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108237", "ses-20221031103712", "pet", ...
                "sub-108237_ses-20221031103712_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(9) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108237", "ses-20221031110638", "pet", ...
                "sub-108237_ses-20221031110638_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(10) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108237", "ses-20221031113804", "pet", ...
                "sub-108237_ses-20221031113804_trc-fdg_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            
            console_clock(11:15) = seconds(26);
            bids_fqfn(11) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108254", "ses-20221116095143", "pet", ...
                "sub-108254_ses-20221116095143_trc-co_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(12) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108254", "ses-20221116100858", "pet", ...
                "sub-108254_ses-20221116100858_trc-oo_proc-delay30-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(13) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108254", "ses-20221116102328", "pet", ...
                "sub-108254_ses-20221116102328_trc-oo_proc-delay30-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(14) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108254", "ses-20221116104751", "pet", ...
                "sub-108254_ses-20221116104751_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(15) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108254", "ses-20221116115244", "pet", ...
                "sub-108254_ses-20221116115244_trc-fdg_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            
            console_clock(16:20) = seconds(43);
            bids_fqfn(16) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108250", "ses-20221207093856", "pet", ...
                "sub-108250_ses-20221207093856_trc-co_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(17) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108250", "ses-20221207095507", "pet", ...
                "sub-108250_ses-20221207095507_trc-oo_proc-delay30-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(18) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108250", "ses-20221207100946", "pet", ...
                "sub-108250_ses-20221207100946_trc-oo_proc-delay30-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(19) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108250", "ses-20221207102944", "pet", ...
                "sub-108250_ses-20221207102944_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(20) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108250", "ses-20221207104909", "pet", ...
                "sub-108250_ses-20221207104909_trc-fdg_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            
            console_clock(21:25) = seconds(25);
            bids_fqfn(21) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108284", "ses-20230220093702", "pet", ...
                "sub-108284_ses-20230220093702_trc-co_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(22) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108284", "ses-20230220095210", "pet", ...
                "sub-108284_ses-20230220095210_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(23) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108284", "ses-20230220101103", "pet", ...
                "sub-108284_ses-20230220101103_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(24) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108284", "ses-20230220103226", "pet", ...
                "sub-108284_ses-20230220103226_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(25) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108284", "ses-20230220112328", "pet", ...
                "sub-108284_ses-20230220112328_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");

            console_clock(26:30) = seconds(32);
            bids_fqfn(26) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108306", "ses-20230227103048", "pet", ...
                "sub-108306_ses-20230227103048_trc-co_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(27) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108306", "ses-20230227104631", "pet", ...
                "sub-108306_ses-20230227104631_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(28) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108306", "ses-20230227112148", "pet", ...
                "sub-108306_ses-20230227112148_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(29) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108306", "ses-20230227113853", "pet", ...
                "sub-108306_ses-20230227113853_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(30) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108306", "ses-20230227115809", "pet", ...
                "sub-108306_ses-20230227115809_trc-fdg_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz"); 

            for row = 1:length(bids_fqfn)
                fp = mybasename(bids_fqfn(row));
                re = regexp(fp, "(?<sub>sub-\d{6})_(?<ses>ses-\d{14})_trc\-(?<tr>[a-z]+)_proc\S+", "names");
                sub(row) = re.sub;
                ses(row) = re.ses;
                tracer(row) = re.tr; %#ok<*AGROW>
            end
            bids_fqfn = ascol(bids_fqfn);
            sub = ascol(sub);
            ses = ascol(ses);       
            tracer = ascol(tracer);

            this.table_maframes_nii_ = table(sub, ses, bids_fqfn, tracer, console_clock);
            T = this.table_maframes_nii_;

            %% local filesystems

            if contains(hostname, "twistor")
                this.strrep_bids_fqfn();
            end
        end
        function T = table_idif_nest(this, matched)
            arguments
                this mlvg.Lee2024
                matched logical = true
            end

            if ~isempty(this.table_idif_nest_)
                T = this.table_idif_nest_;
                return
            end

            T = this.table_maframes_nii;
            T.bids_fqfn = strrep(T.bids_fqfn, "delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz", "MipIdif_idif_dynesty-Boxcar-ideal.csv");
            T.bids_fqfn = strrep(T.bids_fqfn, "sourcedata", "derivatives");

            % include only files on the filesystem
            if matched
                [T0,found] = this.table_idif_nii(true);
                T = T(found, :);
                T = this.add_csv(T);

                T = addvars(T, T0.tidx, NewVariableNames={'tidx'});
                T = this.apply_adjusted_timesMid_img(T, "idif");
                this.table_idif_nest_ = T;
            end
        end
        function [T,found] = table_idif_nii(this, matched)
            arguments
                this mlvg.Lee2024
                matched logical = true
            end

            if ~isempty(this.table_idif_nii_) && ~isempty(this.found_)
                T = this.table_idif_nii_;
                found = this.found_;                
                return
            end

            T = this.table_maframes_nii;
            T.bids_fqfn = strrep(T.bids_fqfn, "delay0-BrainMoCo2-createNiftiMovingAvgFrames", "MipIdif_idif");
            T.bids_fqfn = strrep(T.bids_fqfn, "sourcedata", "derivatives");

            % include only files on the filesystem
            if matched
                [~,found] = this.table_twilite_nii(true);
                T = T(found, :);
                T = this.add_imaging(T);

                T = this.find_adjusted_timesMid_img(T, "idif");
                T = this.apply_adjusted_timesMid_img(T, "idif");
                this.table_idif_nii_ = T;
            end

        end
        function T = table_twilite_nest(this, matched)
            arguments
                this mlvg.Lee2024
                matched logical = true
            end

            if ~isempty(this.table_twilite_nest_)
                T = this.table_twilite_nest_;
                return
            end

            T = this.table_maframes_nii;
            T = this.truncate_ses(T);
            T.bids_fqfn = strrep(T.bids_fqfn, "delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz", "TwiliteKit-do-make-input-func-nomodel_inputfunc_dynesty-RadialArtery-ideal.csv");
            T.bids_fqfn = strrep(T.bids_fqfn, "derivatives", "sourcedata");

            % include only files on the filesystem
            if matched
                [T0,found] = this.table_twilite_nii(true);
                T = T(found, :);
                T = this.add_csv(T);
                
                T = addvars(T, T0.tidx, NewVariableNames={'tidx'});
                T = this.apply_adjusted_timesMid_img(T, "twil");
                this.table_twilite_nest_ = T;
            end
        end
        function [T,found] = table_twilite_nii(this, matched)
            arguments
                this mlvg.Lee2024
                matched logical = true
            end

            if ~isempty(this.table_twilite_nii_) && ~isempty(this.found_)
                T = this.table_twilite_nii_;
                found = this.found_;
                return
            end
            
            T = this.table_maframes_nii;
            T = this.truncate_ses(T);
            T.bids_fqfn = strrep(T.bids_fqfn, "derivatives", "sourcedata");
            T.bids_fqfn = strrep(T.bids_fqfn, "delay0-BrainMoCo2-createNiftiMovingAvgFrames", "TwiliteKit-do-make-input-func-nomodel_inputfunc");
            found = [];

            % include only files on the filesystem
            if matched
                this.found_ = isfile(T.bids_fqfn);
                found = this.found_;
                T = T(found, :);
                T = this.add_imaging(T);     

                T = this.find_adjusted_timesMid_img(T, "twil");
                T = this.apply_adjusted_timesMid_img(T, "twil");
                this.table_twilite_nii_ = T;
            end
        end

        function this = Lee2024(varargin)
        end
    end

    methods (Static)
        function ic = create_kernel(hct, Nt, opts)
            arguments
                hct double
                Nt double = 121
                opts.do_save logical = false
                opts.do_plot logical = false
            end

            cath = mlswisstrace.Catheter_DT20190930(hct=hct);            

            fqfn = fullfile(getenv("HOME"), "Singularity", ...
               "CCIR_01211", "sourcedata", "sub-108293", "ses-20210421155709", "pet", ...
               "sub-108293_ses-20210421155709_trc-fdg_proc-delay0-BrainMoCo2-createNiftiStatic.nii.gz");
            ifc = mlfourd.ImagingFormatContext2(fqfn);
            ifc.img = cath.kernel(Nt=Nt);
            ifc.fqfp = fullfile(getenv("HOME"), "Singularity", ...
               "CCIR_01211", "sourcedata", "kernel_hct="+hct);
            ifc.json_metadata.times = 0:Nt-1;
            ifc.json_metadata.timesMid = 0.5:Nt-0.5;
            ifc.json_metadata.taus = ones(1, Nt);
            ic = mlfourd.ImagingContext2(ifc);   
            if opts.do_save
                save(ic);
            end
            if opts.do_plot
                plot(ic);
            end
        end

        % Record ID,Event Name,Date of session:,Hematocrit
        % 108237,White Matter Hyperintensities (Arm 4: WMH),10/31/2022,43.9
        % 108250,Measuring Aerobic Glycolysis  (Arm 3: MAG),12/7/2022,42.8
        % 108254,Measuring Aerobic Glycolysis  (Arm 3: MAG),11/16/2022,37.9
        % 108283,White Matter Hyperintensities (Arm 4: WMH),2/15/2023,42.1
        % 108284,White Matter Hyperintensities (Arm 4: WMH),2/20/2023,39.7
        % 108293,Measuring Aerobic Glycolysis  (Arm 3: MAG),4/21/2021,46.8
        % 108306,White Matter Hyperintensities (Arm 4: WMH),2/27/2023,41.1

        function rename_nii(fp, fp1, trc)
            fp = myfileprefix(fp);
            re = regexp(fp1, "(?<ss>sub-\d{6}_ses-\d{14})\S*","names");
            fp1 = re.ss+"_trc-"+trc+"_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc";
            ext = [".nii.gz", ".json"];
            for e = ext
                copyfile(fp+e(1), fp1+e(1));
            end
        end
    end
        
    %% PRIVATE

    properties (Access=private)
        found_
        table_centerlines_
        table_idif_nest_
        table_idif_nii_
        table_maframes_nii_
        table_twilite_nest_
        table_twilite_nii_
    end

    methods (Access = private)
        function T = apply_adjusted_timesMid_img(this, T, if_type, opts)
            arguments
                this mlvg.Lee2024
                T table
                if_type {mustBeTextScalar}
                opts.Dtidx double = -15 % shift tidx to accomodate delay of twilite measurements
            end

            assert(any(contains(T.Properties.VariableNames, "tidx")))

            for row = 1:size(T, 1)
                timesMid = T.timesMid{row};
                img = T.img{row};
                tidx = T.tidx(row);

                switch char(if_type)
                    case {'mipidif', 'idif'}
                        img = img(tidx:end);
                        timesMid = timesMid(tidx:end) - timesMid(tidx);
                    case {'twilite', 'twil'}
                        img = img(tidx+opts.Dtidx:end);
                        timesMid = timesMid(tidx+opts.Dtidx:end) - timesMid(tidx+opts.Dtidx);
                    otherwise
                        error("mlvg:RuntimeError", stackstr())
                end

                T.timesMid{row} = timesMid;
                T.img{row} = img;
            end
        end
        function T = find_adjusted_timesMid_img(this, T, if_type)
            arguments
                this mlvg.Lee2024
                T table
                if_type {mustBeTextScalar}
            end

            tidx = nan(size(T.tracer));
            for row = 1:size(T, 1)
                tracer = T.tracer{row};
                img = T.img{row};

                switch char(if_type)
                    case {'mipidif', 'idif'}
                        tidx(row) = mlvg.Lee2024.find_takeoff(img, if_type, tracer); % rising beyond baseline
                    case {'twilite', 'twil'}
                        tidx(row) = mlvg.Lee2024.find_takeoff(img, if_type, tracer); % rising beyond baseline
                    otherwise
                        error("mlvg:RuntimeError", stackstr())
                end
            end

            if ~contains(T.Properties.VariableNames, "tidx")
                T = addvars(T, tidx, NewVariableNames={'tidx'});
            else
                T.tidx = tidx;
            end
        end
        function T = strrep_bids_fqfn(this, opts)
            %% updates this.table_maframes_nii_, replacing arbitrary strings with new strings

            arguments
                this
                opts.s1 {mustBeTextScalar} = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity")
                opts.s2 {mustBeTextScalar} = fullfile(getenv("SINGULARITY_HOME"))
            end
            
            T = this.table_maframes_nii;
            T.bids_fqfn = strrep(T.bids_fqfn, opts.s1, opts.s2);
            this.table_maframes_nii_ = T;
        end
    end

    methods (Static) %, Access = private)
        function T = add_csv(T)
            timesMid = cell(size(T, 1), 1);
            img = cell(size(T, 1), 1);
            for row = 1:size(T, 1)
                csv = readtable(T.bids_fqfn(row));
                timesMid{row} = double(asrow(csv.times));
                img{row} = double(asrow(csv.ideal))/1e3; % kBq/mL
            end
            T = addvars(T, timesMid, img, NewVariableNames={'timesMid', 'img'});
        end
        function T = add_imaging(T)
            imaging = cell(size(T, 1), 1);
            timesMid = cell(size(T, 1), 1);
            img = cell(size(T, 1), 1);
            for row = 1:size(T, 1)
                ic = mlfourd.ImagingContext2(T.bids_fqfn(row));
                ic.selectImagingTool();
                imaging{row} = ic;
                timesMid_ = asrow(ic.json_metadata.timesMid);
                viable = ~isnan(timesMid_);
                timesMid{row} = double(timesMid_(viable));
                img_ = asrow(ic.imagingFormat.img);
                img{row} = double(img_(viable))/1e3; % kBq/mL
            end
            T = addvars(T, imaging, timesMid, img, NewVariableNames={'imaging', 'timesMid', 'img'});
        end
        function tidx = find_takeoff(img, if_type, trc)
            arguments
                img double
                if_type {mustBeTextScalar} % mipidif, twilite, ...
                trc {mustBeTextScalar} 
            end

            switch char(if_type)
                case {'mipidif', 'idif'}
                    thresh = 1e-3;
                case {'twilite', 'twil'}
                    img = img - mean(img(1:10)); % sub baseline
                    switch trc
                        case {'co', 'oc'}
                            thresh = 0.25;
                        case 'oo'
                            thresh = 0.04;
                        case 'ho'
                            thresh = 0.04;
                        case 'fdg'
                            thresh = 0.07; % by inspection
                        otherwise
                            error("mlvg:RuntimeError", stackstr())
                    end
                otherwise
                    error("mlvg:RuntimeError", ...
                        "%s: %s, %s", stackstr(), if_type, trc)
            end
            [~,tidx] = max(img > thresh*max(img));
        end
        function T = truncate_ses(T)
            arguments
                T table
            end

            for bidx = 1:length(T.bids_fqfn)
                fqfn = T.bids_fqfn(bidx);
                parts = split(fqfn, filesep);
                selected = contains(parts, "ses-");
                [~,sidx] = max(selected);
                selected(sidx+1:end) = false;
                ses = parts(selected);
                ses = extractBetween(ses, 1, 12);
                parts(selected) = ses;
                T.bids_fqfn(bidx) = fullfile(filesep, parts{:});
            end
        end
        function dt = ses2datetime(T)
            arguments
                T table
            end

            dt = NaT(size(T, 1));
            for bidx = 1:length(T.bids_fqfn)
                fqfn = T.bids_fqfn(bidx);
                parts = split(fqfn, filesep);
                selected = contains(parts, "ses-");
                [~,sidx] = max(selected);
                selected(sidx+1:end) = false;
                ses = parts(selected);
                dt(bidx) = datetime(extractBetween(ses, 5, 18), InputFormat="yyyyMMddHHmmss");
            end
        end
        function lims = xlim(trc)
            switch char(trc)
                case {'co', 'oc'}
                    lims = [0, 240];
                case 'oo'
                    lims = [0, 120];
                case 'ho'
                    lims = [0, 120];
                case 'fdg'
                    lims = [0, 120];
                
                otherwise
                    error("mlvg:RuntimeError", ...
                        "%s: %s", stackstr(), trc)
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
