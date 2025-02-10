classdef Lee2024 < handle & mlsystem.IHandle
    %% Top-level user interface for softwares supporting 
    %  Lee, et al. to be submitted to Physics in Biology & Medicine, 2024.
    %
    %  See also:  mlvg.Lee2024_livescript.mlx
    %
    %  N.B.: all tables with numerical entries include inverse-efficiency
    %  
    %  Created 03-Jan-2024 21:45:43 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2024 John J. Lee.
    
    methods

        %% Bland-Altman

        function [q_twil,q_idif,abscissa,ordinate, coloring] = bland_altman_for_model(this, opts)
            arguments
                this mlvg.Lee2024
                opts.tracer {mustBeTextScalar} = ""
                opts.parname {mustBeTextScalar} = ""
                opts.subname {mustBeTextScalar} = ""
            end

            %% IDIF
            T_idif = this.table_idif_nest();
            T_idif.bids_fqfn = strrep(T_idif.bids_fqfn, "-ideal.csv", "-quantiles.csv");
            fqfns = T_idif.bids_fqfn;

            qm = NaN(length(fqfns), length(this.par_labels__));
            for scan = 1:length(qm)
                t = readtable(fqfns{scan});
                qm(scan, :) = t.qm';
            end
            T_idif = addvars(T_idif, qm, NewVariableNames="qm");
            T_idif = splitvars(T_idif, "qm", NewVariableNames=this.par_names__);

            %% Twilite
            T_twil = this.table_twilite_nest();
            T_twil.bids_fqfn = strrep(T_twil.bids_fqfn, "-ideal.csv", "-quantiles.csv");
            fqfns = T_twil.bids_fqfn;

            qm = NaN(length(fqfns), length(this.par_labels__));
            for scan = 1:length(qm)
                t = readtable(fqfns{scan});
                qm(scan, :) = t.qm';
            end
            T_twil = addvars(T_twil, qm, NewVariableNames="qm");
            T_twil = splitvars(T_twil, "qm", NewVariableNames=this.par_names__);

            %% gather and normalize; resplit
            T = [T_idif; T_twil];
            for p = this.par_names__
                T.(p(1)) = normalize(T.(p(1))); %, "medianiqr");
            end            
            N = size(T, 1)/2;
            T_idif = T(1:N, :);
            T_twil = T(N+1:end, :);

            coloring = T_twil.coloring;

            if ~isemptytext(opts.tracer)
                T_idif = T_idif(strcmp(T_idif.tracer, opts.tracer), :);
                T_twil = T_twil(strcmp(T_twil.tracer, opts.tracer), :);
            end
            if ~isemptytext(opts.parname)
                q_idif = T_idif.(opts.parname);
                q_twil = T_twil.(opts.parname);
                abscissa = (q_idif + q_twil) / 2;
                ordinate = (q_idif - q_twil); % avoid ./ abscissa for z-scores
                return
            end
            if ~isemptytext(opts.subname)
                T_idif = T_idif(strcmp(T_idif.sub, opts.subname), :);
                T_twil = T_twil(strcmp(T_twil.sub, opts.subname), :);
            end

            %% conveniences
            p_idif = T_idif{:, 13:21}; % [13,14,15,16,17,18,19,20,21]}; % exclude t_0, \tau_2, \tau_3, A, sigma
            p_twil = T_twil{:, 13:21}; % [13,14,15,16,17,18,19,20,21]};

            q_idif = p_idif(:);
            q_twil = p_twil(:);

            abscissa = (q_idif + q_twil) / 2;
            ordinate = (q_idif - q_twil); % avoid ./ abscissa for z-scores
        end

        function plot_bland_altman_pars(this)
            figure; hold on
            colororder(plasma(9))
            for pidx = [3,4,5,6,7,8,9,10,11]
                [~,~,a,o] = this.bland_altman_for_model(parname="p"+pidx);
                s = scatter(a,o, "filled");
                s.SizeData = 100;
                title("Bland-Altman residuals for normalized parameters");
            end
            fontsize(scale=1.5);
            legend(this.par_labels__(4:12));
            hold off

            % figure; hold on
            % colororder(viridis(11))
            % for pidx = [1,2,3,4,5,6,7,8,9,10,11]
            %     [p,q] = this.bland_altman_for_model(parname="p"+pidx);
            %     s = scatter(p,q, "filled");
            %     s.SizeData = 100;
            %     title("AIF \rightarrow IDIF");
            % end
            % fontsize(scale=1.5);
            % legend(this.par_labels__(2:12));
            % hold off
        end
        function plot_bland_altman_subs(this)
            figure; hold on
            colororder(plasma(6))
            for sidx = 1:6
                [~,~,a,o] = this.bland_altman_for_model(subname=this.sub_names__(sidx));
                s = scatter(a,o, "filled");
                s.SizeData = 100;
                title("Bland-Altman residuals for normalized parameters");
            end
            fontsize(scale=1.5);
            legend("participant " + (1:6));
            hold off
        end
        function plot_bland_altman_tracers(this)
            figure; hold on
            colororder(plasma(4))
            for tidx = 1:4
                [~,~,a,o] = this.bland_altman_for_model(tracer=this.tracers__(tidx));
                s = scatter(a,o, "filled");
                s.SizeData = 100;
                title("Bland-Altman residuals for normalized parameters");
            end
            fontsize(scale=1.5);
            legend(this.tracers__);
            hold off
        end
        
        %% centerlines

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
        
        %% tables

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
        function T = table_idif_nest(this, matched)
            arguments
                this mlvg.Lee2024
                matched logical = true
            end

            if ~isempty(this.table_idif_nest_)
                T = this.table_idif_nest_;
                return
            end

            T = this.table_maframes_timeAppend_nii;
            T.bids_fqfn = strrep(T.bids_fqfn, "delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz", "MipIdif_idif_dynesty-Boxcar-ideal.csv");
            T.bids_fqfn = strrep(T.bids_fqfn, "delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-4.nii.gz", "MipIdif_idif_dynesty-Boxcar-ideal.csv");
            T.bids_fqfn = strrep(T.bids_fqfn, "delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-80.nii.gz", "MipIdif_idif_dynesty-Boxcar-ideal.csv");
            T.bids_fqfn = strrep(T.bids_fqfn, "delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-165.nii.gz", "MipIdif_idif_dynesty-Boxcar-ideal.csv");
            T.bids_fqfn = strrep(T.bids_fqfn, "sourcedata", "derivatives");

            % include only files on the filesystem
            if matched
                [~,found] = this.table_idif_nii(true);
                T = T(found, :);
                T = this.add_csv(T, false);
                T = this.add_peak_position(T);

                %T = addvars(T, T0.tidx, NewVariableNames={'tidx'});
                %T = this.apply_adjusted_timesMid_img(T, "idif");
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
            T.bids_fqfn = strrep(T.bids_fqfn, "delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-4", "MipIdif_idif");
            T.bids_fqfn = strrep(T.bids_fqfn, "delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-80", "MipIdif_idif");
            T.bids_fqfn = strrep(T.bids_fqfn, "delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-165", "MipIdif_idif");
            T.bids_fqfn = strrep(T.bids_fqfn, "sourcedata", "derivatives");

            % include only files on the filesystem
            if matched
                [~,found] = this.table_twilite_nii(true);
                T = T(found, :);
            end
            T = this.add_imaging(T, false);

            %T = this.find_adjusted_timesMid_img(T, "idif");
            %T = this.apply_adjusted_timesMid_img(T, "idif");
            this.table_idif_nii_ = T;
        end
        function T = table_maframes_nii(this)
            %% table of moving-average frames (dynamic)

            if ~isempty(this.table_maframes_nii_)
                T = this.table_maframes_nii_;
                return
            end

            hct = NaN(30, 1);
            inveff = NaN(30, 1);
            colors = viridis(6); % pick 6 interpolated colors from viridis colormap

            %console_clock(1:5) = seconds(31);
            hct(1:5) = 46.8;
            inveff(1:5) = 1;
            coloring(1:5) = repmat({colors(1,:)}, [5,1]);
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

            %console_clock(6:10) = seconds(17);
            hct(6:10) = 43.9;
            inveff(6:10) = 1;
            coloring(6:10) = repmat({colors(2,:)}, [5,1]);
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
            
            %console_clock(11:15) = seconds(26);
            hct(11:15) = 37.9;
            inveff(11:15) = 0.9420/1.4967;
            coloring(11:15) = repmat({colors(3,:)}, [5,1]);
            bids_fqfn(11) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108254", "ses-20221116095143", "pet", ...
                "sub-108254_ses-20221116095143_trc-co_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(12) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108254", "ses-20221116100858", "pet", ...
                "sub-108254_ses-20221116100858_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(13) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108254", "ses-20221116102328", "pet", ...
                "sub-108254_ses-20221116102328_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(14) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108254", "ses-20221116104751", "pet", ...
                "sub-108254_ses-20221116104751_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(15) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108254", "ses-20221116115244", "pet", ...
                "sub-108254_ses-20221116115244_trc-fdg_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            
            %console_clock(16:20) = seconds(43);
            hct(16:20) = 42.8;
            inveff(16:20) = 0.9420/1.6809;
            coloring(16:20) = repmat({colors(4,:)}, [5,1]);
            bids_fqfn(16) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108250", "ses-20221207093856", "pet", ...
                "sub-108250_ses-20221207093856_trc-co_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(17) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108250", "ses-20221207095507", "pet", ...
                "sub-108250_ses-20221207095507_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(18) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108250", "ses-20221207100946", "pet", ...
                "sub-108250_ses-20221207100946_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(19) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108250", "ses-20221207102944", "pet", ...
                "sub-108250_ses-20221207102944_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(20) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108250", "ses-20221207104909", "pet", ...
                "sub-108250_ses-20221207104909_trc-fdg_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            
            %console_clock(21:25) = seconds(25);
            hct(21:25) = 39.7;
            inveff(21:25) = 0.9420/1.7463;
            coloring(21:25) = repmat({colors(5,:)}, [5,1]);
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

            %console_clock(26:30) = seconds(32);
            hct(26:30) = 41.1;
            inveff(26:30) = 0.9420/1.7187;
            coloring(26:30) = repmat({colors(6,:)}, [5,1]);
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

            %console_clock = ascol(console_clock);
            hct = ascol(hct);

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
            hct = ascol(hct);  
            %console_clock = ascol(console_clock);  
            inveff = ascol(inveff);  
            coloring = ascol(coloring);

            this.table_maframes_nii_ = table(sub, ses, bids_fqfn, tracer, hct, inveff, coloring); % , console_clock
            T = this.table_maframes_nii_;

            %% local filesystems

            if contains(hostname, "twistor")
                this.table_maframes_nii_ = this.strrep_bids_fqfn(this.table_maframes_nii_);
                T = this.strrep_bids_fqfn(T);
            end
        end
        function T = table_maframes_timeAppend_nii(this)
            %% table of moving-average frames (dynamic)

            if ~isempty(this.table_maframes_timeAppend_nii_)
                T = this.table_maframes_timeAppend_nii_;
                return
            end

            hct = NaN(30, 1);
            inveff = NaN(30, 1);
            colors = viridis(6); % pick 6 interpolated colors from viridis colormap

            %console_clock(1:5) = seconds(31);
            hct(1:5) = 46.8;
            inveff(1:5) = 1;
            coloring(1:5) = repmat({colors(1,:)}, [5,1]);
            bids_fqfn(1) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108293", "ses-20210421144815", "pet", ...
                "sub-108293_ses-20210421144815_trc-co_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(2) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108293", "ses-20210421150523", "pet", ...
                "sub-108293_ses-20210421150523_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-4.nii.gz");
            bids_fqfn(3) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108293", "ses-20210421152358", "pet", ...
                "sub-108293_ses-20210421152358_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(4) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108293", "ses-20210421154248", "pet", ...
                "sub-108293_ses-20210421154248_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-4.nii.gz");
            bids_fqfn(5) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
               "CCIR_01211", "sourcedata", "sub-108293", "ses-20210421155709", "pet", ...
               "sub-108293_ses-20210421155709_trc-fdg_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-165.nii.gz");           

            %console_clock(6:10) = seconds(17);
            hct(6:10) = 43.9;
            inveff(6:10) = 1;
            coloring(6:10) = repmat({colors(2,:)}, [5,1]);
            bids_fqfn(6) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108237", "ses-20221031100910", "pet", ...
                "sub-108237_ses-20221031100910_trc-co_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(7) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108237", "ses-20221031102320", "pet", ...
                "sub-108237_ses-20221031102320_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-4.nii.gz");
            bids_fqfn(8) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108237", "ses-20221031103712", "pet", ...
                "sub-108237_ses-20221031103712_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-4.nii.gz");
            bids_fqfn(9) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108237", "ses-20221031110638", "pet", ...
                "sub-108237_ses-20221031110638_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(10) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108237", "ses-20221031113804", "pet", ...
                "sub-108237_ses-20221031113804_trc-fdg_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-165.nii.gz");
            
            %console_clock(11:15) = seconds(26);
            hct(11:15) = 37.9;
            inveff(11:15) = 0.9420/1.4967;
            coloring(11:15) = repmat({colors(3,:)}, [5,1]);
            bids_fqfn(11) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108254", "ses-20221116095143", "pet", ...
                "sub-108254_ses-20221116095143_trc-co_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(12) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108254", "ses-20221116100858", "pet", ...
                "sub-108254_ses-20221116100858_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-4.nii.gz");
            bids_fqfn(13) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108254", "ses-20221116102328", "pet", ...
                "sub-108254_ses-20221116102328_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-4.nii.gz");
            bids_fqfn(14) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108254", "ses-20221116104751", "pet", ...
                "sub-108254_ses-20221116104751_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(15) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108254", "ses-20221116115244", "pet", ...
                "sub-108254_ses-20221116115244_trc-fdg_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-165.nii.gz");
            
            %console_clock(16:20) = seconds(43);
            hct(16:20) = 42.8;
            inveff(16:20) = 0.9420/1.6809;
            coloring(16:20) = repmat({colors(4,:)}, [5,1]);
            bids_fqfn(16) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108250", "ses-20221207093856", "pet", ...
                "sub-108250_ses-20221207093856_trc-co_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(17) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108250", "ses-20221207095507", "pet", ...
                "sub-108250_ses-20221207095507_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-4.nii.gz");
            bids_fqfn(18) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108250", "ses-20221207100946", "pet", ...
                "sub-108250_ses-20221207100946_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-4.nii.gz");
            bids_fqfn(19) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108250", "ses-20221207102944", "pet", ...
                "sub-108250_ses-20221207102944_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(20) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108250", "ses-20221207104909", "pet", ...
                "sub-108250_ses-20221207104909_trc-fdg_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-165.nii.gz");
            
            %console_clock(21:25) = seconds(25);
            hct(21:25) = 39.7;
            inveff(21:25) = 0.9420/1.7463;
            coloring(21:25) = repmat({colors(5,:)}, [5,1]);
            bids_fqfn(21) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108284", "ses-20230220093702", "pet", ...
                "sub-108284_ses-20230220093702_trc-co_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(22) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108284", "ses-20230220095210", "pet", ...
                "sub-108284_ses-20230220095210_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-4.nii.gz");
            bids_fqfn(23) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108284", "ses-20230220101103", "pet", ...
                "sub-108284_ses-20230220101103_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-4.nii.gz");
            bids_fqfn(24) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108284", "ses-20230220103226", "pet", ...
                "sub-108284_ses-20230220103226_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(25) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108284", "ses-20230220112328", "pet", ...
                "sub-108284_ses-20230220112328_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");

            %console_clock(26:30) = seconds(32);
            hct(26:30) = 41.1;
            inveff(26:30) = 0.9420/1.7187;
            coloring(26:30) = repmat({colors(6,:)}, [5,1]);
            bids_fqfn(26) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108306", "ses-20230227103048", "pet", ...
                "sub-108306_ses-20230227103048_trc-co_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(27) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108306", "ses-20230227104631", "pet", ...
                "sub-108306_ses-20230227104631_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-80.nii.gz");
            bids_fqfn(28) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108306", "ses-20230227112148", "pet", ...
                "sub-108306_ses-20230227112148_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-4.nii.gz");
            bids_fqfn(29) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108306", "ses-20230227113853", "pet", ...
                "sub-108306_ses-20230227113853_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");
            bids_fqfn(30) = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", ...
                "CCIR_01211", "sourcedata", "sub-108306", "ses-20230227115809", "pet", ...
                "sub-108306_ses-20230227115809_trc-fdg_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-165.nii.gz"); 

            %console_clock = ascol(console_clock);
            hct = ascol(hct);

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
            hct = ascol(hct);  
            %console_clock = ascol(console_clock);  
            inveff = ascol(inveff);  
            coloring = ascol(coloring);

            this.table_maframes_timeAppend_nii_ = table(sub, ses, bids_fqfn, tracer, hct, inveff, coloring); % , console_clock
            T = this.table_maframes_timeAppend_nii_;

            %% local filesystems

            if contains(hostname, "twistor")
                this.table_maframes_timeAppend_nii_ = this.strrep_bids_fqfn(this.table_maframes_timeAppend_nii_);
                T = this.strrep_bids_fqfn(T);
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

            T = this.table_maframes_timeAppend_nii;
            T.bids_fqfn = strrep(T.bids_fqfn, "delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz", "TwiliteKit-do-make-input-func-nomodel_inputfunc_dynesty-RadialArtery-ideal.csv");
            T.bids_fqfn = strrep(T.bids_fqfn, "delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-4.nii.gz", "TwiliteKit-do-make-input-func-nomodel_inputfunc_dynesty-RadialArtery-ideal.csv");
            T.bids_fqfn = strrep(T.bids_fqfn, "delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-80.nii.gz", "TwiliteKit-do-make-input-func-nomodel_inputfunc_dynesty-RadialArtery-ideal.csv");
            T.bids_fqfn = strrep(T.bids_fqfn, "delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-165.nii.gz", "TwiliteKit-do-make-input-func-nomodel_inputfunc_dynesty-RadialArtery-ideal.csv");
            T.bids_fqfn = strrep(T.bids_fqfn, "sourcedata", "derivatives");

            % include only files on the filesystem
            if matched
                [~,found] = this.table_twilite_nii(true);
                T = T(found, :);
                T = this.add_csv(T, true);
                T = this.add_peak_position(T);
                
                %T = addvars(T, T0.tidx, NewVariableNames={'tidx'});
                %T = this.apply_adjusted_timesMid_img(T, "twil");
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
            T.bids_fqfn = strrep(T.bids_fqfn, "delay0-BrainMoCo2-createNiftiMovingAvgFrames", "TwiliteKit-do-make-input-func-nomodel_inputfunc");
            T.bids_fqfn = strrep(T.bids_fqfn, "delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-4", "TwiliteKit-do-make-input-func-nomodel_inputfunc");
            T.bids_fqfn = strrep(T.bids_fqfn, "delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-80", "TwiliteKit-do-make-input-func-nomodel_inputfunc");
            T.bids_fqfn = strrep(T.bids_fqfn, "delay0-BrainMoCo2-createNiftiMovingAvgFrames_timeAppend-165", "TwiliteKit-do-make-input-func-nomodel_inputfunc");
            T.bids_fqfn = strrep(T.bids_fqfn, "sourcedata", "derivatives");
            found = [];

            % include only files on the filesystem
            if matched
                this.found_ = isfile(T.bids_fqfn);
                found = this.found_;
                T = T(found, :);
                T = this.add_imaging(T, true);

                %T = this.find_adjusted_timesMid_img(T, "twil");
                %T = this.apply_adjusted_timesMid_img(T, "twil");
                this.table_twilite_nii_ = T;
            end
        end

        %% bigraphs

        function build_all_twilite(this)
            T = this.table_maframes_nii();
            for row = 17:18 % 1:size(T, 1)
                try
                    fqfn = T.bids_fqfn{row};
                    [pth,fp] = myfileparts(fqfn);
                    TRC = this.filename2tracer(fqfn);

                    % registry
                    reg = mlvg.Ccir1211Registry.instance();
                    reg.tracer = upper(TRC);

                    % crv
                    crv_fqfn = fullfile(getenv("HOME"), ...
                        "Documents", "CCIRRadMeasurements", "Twilite", "CRV", this.fqfn2crv(fqfn));

                    % ifk
                    ifk = mlkinetics.InputFuncKit.create_from_tags( ...
                        bids_fqfn=fqfn, ...
                        bids_tags="ccir1211", ...
                        tracer_tags=TRC, ...
                        scanner_tags="vision", ...
                        input_func_tags="twilite-nomodel-nocrop", ...
                        input_func_fqfn=crv_fqfn, ...
                        hct=T.hct(row)); % '-nocrop' ~ include Twilite activity prior to start of scanner
                    twil_ic = ifk.do_make_activity_density(doDecayCorrection=true);
                    plot(twil_ic);
                    twil_ic.filepath = strrep(pth, "sourcedata", "derivatives");
                    twil_ic.fileprefix = strrep(fp, ...
                        "delay0-BrainMoCo2-createNiftiMovingAvgFrames", ...
                        "TwiliteKit-do-make-input-func-nomodel_inputfunc");
                    twil_ic.fileprefix = strrep(fp, ...
                        "delay30-BrainMoCo2-createNiftiMovingAvgFrames", ...
                        "TwiliteKit-do-make-input-func-nomodel_inputfunc");
                    twil_ic.save();
                catch ME
                    handwarning(ME)
                end
            end
        end
        function deconv_twilite(this)
            T = this.table_twilite_nii();
            for row = 1:size(T, 1)
                fqfn = T.bids_fqfn{row};
                fqfp = myfileprefix(fqfn);
                ic = mlfourd.ImagingFormatContext2(fqfn);
                t = ic.json_metadata.timesMid;
                rho = ic.img;
                plot(t, rho); hold on

                cath = mlswisstrace.Catheter_DT20190930( ...
                    Measurement=rho, ...
                    hct=T.hct(row), ...
                    tracer=T.tracer{row}, ...
                    fqfileprefix=fqfp+"_deconv_twilite");
                rho1 = cath.deconv();
                plot(t, rho1); hold off

                title(ic.fqfn)
            end
        end
        
        function call_ifk(this, row, opts)
            arguments
                this mlvg.Lee2024
                row double
                opts.steps logical = [0,0,0,0,1,0]
                opts.delete_large_files logical = false
                opts.reference_tracer {mustBeTextScalar} = "fdg"
                opts.frac_select double = mlaif.MipIdif.ALPHA
                opts.dilate_m double = 0
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
                    steps=opts.steps, ...
                    reference_tracer=opts.reference_tracer, ...
                    delete_large_files=opts.delete_large_files, ...
                    frac_select=opts.frac_select, ...
                    dilate_m=opts.dilate_m);
            catch ME
                handwarning(ME)
            end
        end
        
        function plot_bigraph(this, trc, ti)
            Ttwil = this.table_twilite_nii();
            Tidif = this.table_idif_nii();

            Utwil = Ttwil(Ttwil.tracer == trc, :)
            Uidif = Tidif(Tidif.tracer == trc, :)

            figure; hold on
            for r = 1:size(Uidif, 1)
                timesMid__ = Uidif.timesMid{r};
                img__ = Uidif.img{r};
                if strcmp(trc, "oo") && r == 6
                    selected = Uidif.img{r} > 2;
                    selected(1:10) = true;
                    timesMid__ = timesMid__(selected);
                    img__ = img__(selected);
                end
                plot(timesMid__, img__, Color=Uidif.coloring{r}, LineWidth=4);
            end
            hold off
            set(gca, XTickLabel=[]);
            if strcmp(trc, "co")
                ylabel("activity (kBq/mL)", FontSize=18); end
            fontsize(scale = 2)
            %title(ti, FontSize=24)
            annotation('textbox', [.55 .58 .3 .3], 'String', ti, 'FitBoxToText', 'on', 'LineStyle', 'none', 'FontSize', 42); 
            xlim(this.xlim(trc));
            ylim(this.ylim(trc));

            figure; hold on
            for r = 1:size(Utwil, 1)
                plot(Utwil.timesMid{r}, Utwil.img{r}, Color=Utwil.coloring{r}, LineWidth=4); end; hold off
            set(gca, 'YDir', 'reverse')
            xlabel("time (s)", FontSize=18);
            fontsize(scale = 2)
            xlim(this.xlim(trc));
            ylim(this.ylim(trc));
        end
        function plot_bigraph_nest(this, trc, ti)
            Ttwil = this.table_twilite_nest();
            Tidif = this.table_idif_nest();

            Utwil = Ttwil(Ttwil.tracer == trc, :)
            Uidif = Tidif(Tidif.tracer == trc, :)

            figure; hold on
            for r = 1:size(Uidif, 1)
                timesMid__ = Uidif.timesMid{r};
                img__ = Uidif.img{r};
                % if strcmp(trc, "oo") && r == 6
                %     selected = Uidif.img{r} > 2;
                %     selected(1:10) = true;
                %     timesMid__ = timesMid__(selected);
                %     img__ = img__(selected);
                % end
                plot(timesMid__, img__, Color=Uidif.coloring{r}, LineWidth=4); 
            end
            hold off
            set(gca, XTickLabel=[]);
            if strcmp(trc, "co")
                ylabel("activity (kBq/mL)", FontSize=18); end
            fontsize(scale = 2)
            %title(ti, FontSize=24)
            annotation('textbox', [.55 .58 .3 .3], 'String', ti, 'FitBoxToText', 'on', 'LineStyle', 'none', 'FontSize', 42); 
            xlim(this.xlim(trc));
            ylim(this.ylim(trc, true));

            figure; hold on
            for r = 1:size(Utwil, 1)
                plot(Utwil.timesMid{r}, Utwil.img{r}, Color=Utwil.coloring{r}, LineWidth=4); end; hold off
            set(gca, 'YDir', 'reverse')
            xlabel("time (s)", FontSize=18)
            fontsize(scale = 2)
            xlim(this.xlim(trc));
            ylim(this.ylim(trc, true));
        end
        
        %% build recovery coefficients

        function C = build_all_recovery_coeff(this)
            U = this.table_idif_nii();
            U = U(U.tracer == "co", :);

            V = this.table_twilite_nii();
            V = V(V.tracer == "co", :);

            C = cell(1, size(V, 1));
            rc = [];
            idif = [];
            twil = [];
            residual = [];
            estimator = [];
            rgb = [];
            for row = 1:size(U, 1)
                ic = U.imaging{row};
                coloring = U.coloring{row};
                [ifc_rc, ifc_idif, ifc_twil] = this.build_recovery_coeff(U(row,:), V(row,:), ic);
                C{row} = struct( ...
                    "rc_img", ifc_rc.img, ...
                    "idif_img", ifc_idif.img, ...
                    "twil_img", ifc_twil.img, ...
                    "residual", ifc_idif.img - ifc_twil.img, ...
                    "estimator", 0.5*(ifc_idif.img + ifc_twil.img), ...
                    "coloring", coloring);

                rc = [rc, asrow(ifc_rc.img)];
                idif = [idif, asrow(ifc_idif.img)];
                twil = [twil, asrow(ifc_twil.img)];
                residual = [residual, asrow(ifc_idif.img - ifc_twil.img)];
                estimator = [estimator, asrow(0.5*(ifc_idif.img + ifc_twil.img))];
                rgb = [rgb; repmat(coloring, flip(size(ifc_rc)))];
            end

            ONES = ones(size(rc));
            LS = linspace(0, 250, 500);

            figure;
            hold on
            plot(LS, LS, "-", LineWidth=0.25, Color=[0.5, 0.5, 0.5])
            s1 = scatter(twil, idif, 100*ONES, rgb, "filled");
            s1.AlphaData = 0.05*ONES;
            s1.MarkerFaceAlpha = 'flat';
            title("Bland-Altman scatter for steady-state [^{15}O]CO")
            xlabel("AIF activity (kBq/mL)")
            ylabel("IDIF activity (kBq/mL)")
            xlim([0 220])
            ylim([0 220])
            fontsize(scale=1.5)
            hold off

            figure;
            hold on
            plot(LS, zeros(size(LS)), "-", LineWidth=0.25, Color=[0.5, 0.5, 0.5])
            s2 = scatter(estimator, residual, 100*ONES, rgb, "filled");
            s2.AlphaData = 0.05*ONES;
            s2.MarkerFaceAlpha = 'flat';
            title("Bland-Altman residuals for steady-state [^{15}O]CO")   
            xlabel("unbiased estimator of activity (kBq/mL)")
            ylabel("IDIF activity - AIF activity (kBq/mL)")
            xlim([0 160])
            %ylim([0 250])
            fontsize(scale=1.5)
            hold off

            figure;
            hold on
            plot(LS, median(rc)*ones(size(LS)), "-", LineWidth=0.25, Color=[0.5, 0.5, 0.5])
            s3 = scatter(estimator, rc, 100*ONES, rgb, "filled");
            s3.AlphaData = 0.05*ONES;
            s3.MarkerFaceAlpha = 'flat';
            title("Recovery coefficients for steady-state [^{15}O]CO") 
            xlabel("unbiased estimator of activity (kBq/mL)")
            ylabel("recovery coefficient")
            xlim([0 160])
            %ylim([0 250])
            fontsize(scale=1.5)
            hold off

            fprintf("median(rc) -> %g\n", median(rc))

        end
        function C = build_all_recovery_coeff2(this, tracer)
            arguments
                this mlvg.Lee2024
                tracer {mustBeTextScalar}
            end

            U = this.table_idif_nest();
            U = U(U.tracer == tracer, :);
            V = this.table_twilite_nest();
            V = V(V.tracer == tracer, :);

            C = cell(1, size(V, 1));
            rc = [];
            idif = [];
            twil = [];
            residual = [];
            estimator = [];
            rgb = [];
            for row = 1:size(U, 1)
                coloring = U.coloring{row};
                [ifc_rc, ifc_idif, ifc_twil] = this.build_recovery_coeff2(U(row,:), V(row,:));
                C{row} = struct( ...
                    "rc_img", ifc_rc.img, ...
                    "idif_img", ifc_idif.img, ...
                    "twil_img", ifc_twil.img, ...
                    "residual", ifc_idif.img - ifc_twil.img, ...
                    "estimator", 0.5*(ifc_idif.img + ifc_twil.img), ...
                    "coloring", coloring);

                rc = [rc, asrow(ifc_rc.img)];
                idif = [idif, asrow(ifc_idif.img)];
                twil = [twil, asrow(ifc_twil.img)];
                residual = [residual, asrow(ifc_idif.img - ifc_twil.img)];
                estimator = [estimator, asrow(0.5*(ifc_idif.img + ifc_twil.img))];
                rgb = [rgb; repmat(coloring, flip(size(ifc_rc)))];
            end

            ONES = ones(size(rc));

            figure;
            hold on
            LS = linspace(0, max(estimator), 500);
            plot(LS, zeros(size(LS)), "-", LineWidth=0.25, Color=[0.5, 0.5, 0.5])
            s2 = scatter(estimator, residual, 100*ONES, rgb, "filled");
            s2.AlphaData = 0.05*ONES;
            s2.MarkerFaceAlpha = 0.2;
            title("Bland-Altman residuals for " + upper(tracer))
            xlabel("unbiased estimator of decays (kBq s/mL)")
            ylabel("IDIF decays - AIF decays (kBq s/mL)")
            %xlim([0 30e3])
            %ylim([min() max()])
            fontsize(scale=1.5)
            hold off

            return

            figure;
            hold on
            LS = linspace(0, max(twil), 500);
            plot(LS, LS, "-", LineWidth=0.25, Color=[0.5, 0.5, 0.5])
            s1 = scatter(twil, idif, 100*ONES, rgb, "filled");
            s1.AlphaData = 0.05*ONES;
            s1.MarkerFaceAlpha = 0.2;
            title("Bland-Altman scatter for " + upper(tracer))
            xlabel("AIF decays (kBq s/mL)")
            ylabel("IDIF decays (kBq s/mL)")
            %xlim([0 220])
            %ylim([0 220])
            fontsize(scale=1.5)
            hold off

            figure;
            hold on
            LS = linspace(0, max(estimator), 500);
            plot(LS, ones(size(LS)), "-", LineWidth=0.25, Color=[0.5, 0.5, 0.5])
            s3 = scatter(estimator, rc, 100*ONES, rgb, "filled");
            s3.AlphaData = 0.05*ONES;
            s3.MarkerFaceAlpha = 0.2;
            title("Recovery coefficients for " + upper(tracer))
            xlabel("unbiased estimator of decays (kBq s/mL)")
            ylabel("recovery coefficient")
            %xlim([0 160])
            ylim([0 2])
            fontsize(scale=1.5)
            hold off

            fprintf("median(rc) -> %g\n", median(rc))

        end
        function [ifc_rc, ifc_idif, ifc_twil] = build_recovery_coeff2(this, U, V)

            if contains(U.tracer, "fdg")
                [ifc_rc, ifc_idif, ifc_twil] = this.build_recovery_coeff2_fdg(U, V);
                return
            end

            halflife = 122.2416;  % sec
            U_ic = mlfourd.ImagingContext2(myfileprefix(U.bids_fqfn{1}) + ".nii.gz");
            ifc_idif = U_ic.imagingFormat;
            ifc_twil = U_ic.imagingFormat;
            ifc_rc = U_ic.imagingFormat;

            %U_tM = U.timesMid{1};
            %V_tM = V.timesMid{1};
            U_img = U.img{1};
            V_img = V.img{1};
            U_pp = U.peak_position(1);
            V_pp = V.peak_position(1);
            Dt = V_pp - U_pp;
            V_img = V_img .* 2.^(Dt/ halflife);
            Nt = min(length(U_img) - U_pp + 1, length(V_img) - V_pp + 1);
            ifc_idif.img = nan(1, Nt);
            ifc_twil.img = nan(1, Nt);
            ifc_rc.img = nan(1, Nt);

            for t = 1:Nt
                ifc_idif.img(t) = sum(U_img(1:U_pp+t-1));
                ifc_twil.img(t) = sum(V_img(1:V_pp+t-1));
            end
            ifc_rc.img = ifc_idif.img ./ ifc_twil.img;

            s = split(U_ic.fileprefix, "proc-");
            ifc_idif.fileprefix = s(1) + "proc-build-rc2_idif";
            ifc_idif.save();
            ifc_twil.fileprefix = s(1) + "proc-build-rc2_twil";
            ifc_twil.save();
            ifc_rc.fileprefix = s(1) + "proc-build-rc2_rc";
            ifc_rc.save();
        end
        function [ifc_rc, ifc_idif, ifc_twil] = build_recovery_coeff2_fdg(~, U, V)
            %% Specifically manage FDG.  
            %  Args:
            %    U (table): information re. IDIFs.
            %               Measurements must be decay-corrected.
            %    V (table): information re. AIFs
            %               Measurements must be decaying, viz., decay corrections take
            %               counting measurements to time of arterial sampling
            %  Returns:
            %    mlfourd.ImagingFormatContext2: recovery coefficients


            halflife = 1.82951 * 3600;  % sec
            U_ic = mlfourd.ImagingContext2(myfileprefix(U.bids_fqfn{1}) + "-embed.nii.gz");
            V_ic = mlfourd.ImagingContext2(myfileprefix(V.bids_fqfn{1}) + "-embed.nii.gz");
            ifc_idif = U_ic.imagingFormat;
            ifc_twil = V_ic.imagingFormat;
            ifc_rc = copy(U_ic.imagingFormat);

            %U_tM = U.timesMid{1};
            %V_tM = V.timesMid{1};
            U_img = ifc_idif.img;
            V_img = ifc_twil.img;
            [~,U_pp] = max(U_img);
            [~,V_pp] = max(V_img);
            Dt = V_pp - U_pp;
            V_img = V_img .* 2.^(Dt/ halflife);
            Nt = min(length(U_img) - U_pp + 1, length(V_img) - V_pp + 1);
            ifc_idif.img = nan(1, Nt);
            ifc_twil.img = nan(1, Nt);
            ifc_rc.img = nan(1, Nt);

            for t = 1:Nt
                ifc_idif.img(t) = sum(U_img(1:U_pp+t-1));
                ifc_twil.img(t) = sum(V_img(1:V_pp+t-1));
            end
            ifc_rc.img = ifc_idif.img ./ ifc_twil.img;

            s = split(U_ic.fileprefix, "proc-");
            ifc_idif.fileprefix = s(1) + "proc-build-rc2_idif";
            ifc_idif.save();
            ifc_twil.fileprefix = s(1) + "proc-build-rc2_twil";
            ifc_twil.save();
            ifc_rc.fileprefix = s(1) + "proc-build-rc2_rc";
            ifc_rc.save();
        end
        function [ifc_rc, ifc_idif, ifc_twil] = build_recovery_coeff(~, U, V, ic)

            U_img = U.img{1};
            V_img = V.img{1};
            U_tM = U.timesMid{1}; % - U.timesMid{1}(1);
            V_tM = V.timesMid{1}; % - V.timesMid{1}(1);
            %halflife = 122.2416;

            %U_img = U_img .* 2.^(-U_tM/halflife);
            %V_img = V_img .* 2.^(-V_tM/halflife);

            U_select = U_tM > V_tM(1) + 120;
            U_img = U_img(U_select);
            V_select = V_tM > V_tM(1) + 120;
            V_img = V_img(V_select);
            Nt = min(length(U_img), length(V_img));

            ifc_idif = ic.imagingFormat;
            ifc_idif.img = U_img(1:Nt);
            ifc_idif.fileprefix = strrep(ic.fileprefix, "MipIdif_idif", stackstr(use_dashes=true) + "_idif");
            %figure; plot(ifc_idif)
            ifc_idif.save();

            ifc_twil = ic.imagingFormat;
            ifc_twil.img = V_img(1:Nt);
            ifc_twil.fileprefix = strrep(ic.fileprefix, "MipIdif_idif", stackstr(use_dashes=true) + "_twil");
            %figure; plot(ifc_twil)
            ifc_twil.save();

            ifc_rc = ic.imagingFormat;
            ifc_rc.img = U_img(1:Nt) ./ V_img(1:Nt);
            ifc_rc.fileprefix = strrep(ic.fileprefix, "MipIdif_idif", stackstr(use_dashes=true) + "_rc");
            %figure; plot(ifc_rc)
            ifc_rc.save();
        end
        
        %% more building
        
        function build_all_maframes_timeAppend_parc(this)
            import mlkinetics.*

            T = this.table_maframes_timeAppend_nii;
            parfor (row = 1:size(T, 1), 8)
                petFqfn = T.bids_fqfn(row);
                this.build_maframes_timeAppend_parc(petFqfn);
            end
        end
        function ic1 = build_maframes_timeAppend_parc(this, petFqfn, opts)
            arguments
                this mlvg.Lee2024
                petFqfn {mustBeFile}
                opts.tag string = ""
            end
            if ~isemptytext(opts.tag) && ~startsWith(opts.tag, "-")
                opts.tag = "-" + opts.tag;
            end

            import mlkinetics.*

            petMed = mlvg.Ccir1211Mediator.create(petFqfn);

            omat = fullfile(petMed.derivPetPath, "T1w_on_" + petMed.imagingReference.fileprefix + ".mat");
            omat_no_tag = this.remove_tau_tag(omat);
            if ~isfile(omat) && isfile(omat_no_tag)
                omat = omat_no_tag;
            end
            flirt = mlfsl.Flirt( ...
                'in', petMed.schaeffer_ic, ...
                'ref', petMed.imagingReference, ...
                'out', strcat(petMed.fqfp, "-schaeffer.nii.gz"), ...
                'omat', omat, ...
                'bins', 1024, ...
                'interp', 'nearestneighbour', ...
                'noclobber', false);
            flirt.applyXfm();

            schBK = BidsKit.create( ...
                bids_tags="ccir1211", bids_fqfn=flirt.out.fqfn);
            rk = RepresentationKit.create( ...
                representation_tags="trivial");
            pk = ParcKit.create( ...
                bids_kit=schBK, representation_kit=rk, parc_tags="schaeffer-schaeffer");

            p = pk.make_parc();
            ic1 = p.reshape_to_parc(petMed.imagingContext);
            disp(ic1)
            %ic1.view()
            %ic1.view()
            ic1.filepath = strrep(ic1.filepath, "sourcedata", "derivatives");
            ic1.fileprefix = ic1.fileprefix + opts.tag;
            ic1.save();
        end
        function build_plasma_fdg_input_func(this)
            %% applies mlraichle.RBCPartition.wb2plasma() for plasma corrections needed by FDG

            idif_to_glob = "sub-*_ses-*_trc-fdg_proc-MipIdif_idif_dynesty-Boxcar-ideal-embed.nii.gz";
            twil_to_glob = "sub-*_ses-*_trc-fdg_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc_dynesty-RadialArtery-ideal-embed.nii.gz";
            for tg = [idif_to_glob, twil_to_glob]
                mg = mglob(fullfile( ...
                    getenv("SINGULARITY_HOME"), ...
                    "CCIR_01211", ...
                    "derivatives", ...
                    "sub-*", "ses-*", "pet", tg));
                for an_mg = mg
                    ic =  mlfourd.ImagingContext2(an_mg);
                    fp_ = ic.fileprefix;
                    ic = mlraichle.RBCPartition.wb2plasma( ...
                        ic, ...
                        this.found_hct(an_mg));
                    ic.fileprefix = strrep(fp_, "-embed", "-plasma");
                    save(ic);
                end
            end
        end
        function build_qms_as_dtseries(this, filenames_toglob, tags)
            arguments
                this mlvg.Lee2024
                filenames_toglob {mustBeText}  % appended to $SINGULARITY_HOME/CCIR_01211/derivatives/sub-*/ses-*/pet
                tags {mustBeTextScalar}
            end

            for tog = asrow(filenames_toglob)
                globbed = mglob( ...
                    fullfile( ...
                        getenv("SINGULARITY_HOME"), "CCIR_01211", "derivatives", "sub-*", "ses-*", "pet", ...
                        tog(1)));
                trc = this.filename2tracer(globbed(1));
                new_fqfn = fullfile( ...
                    getenv("SINGULARITY_HOME"), "CCIR_01211", "derivatives", ...
                    sprintf("sub-all_ses-all_trc-%s_proc-schaeffer-%s", trc, tags));
                Nkind = 1; % this.filename2Nkind(globbed(1));
                for k = 1:Nkind
                    m = this.infer_multiplier(globbed(1), k);
                    cifti = mlvg.Lee2024.qms_to_dtseries(globbed, new_fqfn+"-kind"+k, qm_kind=k, multiplier=m);

                    %% write median(, 2)

                    cifti.cdata = median(cifti.cdata, 2);
                    cifti.diminfo{2} = cifti_diminfo_make_scalars(1);
                    [pth,fp] = myfileparts(new_fqfn+"-kind"+k);
                    cifti_write(cifti, convertStringsToChars(fullfile(pth, fp + "_median_2.dscalar.nii")));
                end
            end
        end       
        function build_RadialArtery_ideal_recalibrated(this)
            T = this.table_maframes_timeAppend_nii();
            for row = 1:size(T, 1)
                fqfn = T.bids_fqfn(row);
                fqfn = strrep(fqfn, "sourcedata", "derivatives");
                idx = strfind(fqfn, "_proc-");
                fqfnc = convertStringsToChars(fqfn);
                fqfn_twil = fqfnc(1:idx) + "proc-TwiliteKit-do-make-input-func-nomodel_inputfunc_dynesty-RadialArtery-ideal.nii.gz";
                fqfn_recal = fqfnc(1:idx) + "proc-TwiliteKit-do-make-input-func-nomodel_inputfunc_dynesty-RadialArtery-ideal-recalibrated.nii.gz";
                
                if ~isfile(fqfn_twil)
                    continue
                end
                disp(fqfn_twil)
                ic = mlfourd.ImagingContext2(fqfn_twil);
                ic = ic ./ T.inveff(row);
                ic.fqfn = fqfn_recal;
                ic.save()
            end
        end
        function build_TwiliteKit_nomodel_recalibrated_inputfunc(this)
            T = this.table_maframes_timeAppend_nii();
            for row = 1:size(T, 1)
                fqfn = T.bids_fqfn(row);
                fqfn = strrep(fqfn, "sourcedata", "derivatives");
                idx = strfind(fqfn, "_proc-");
                fqfnc = convertStringsToChars(fqfn);
                fqfn_twil = fqfnc(1:idx) + "proc-TwiliteKit-do-make-input-func-nomodel_inputfunc.nii.gz";
                fqfn_recal = fqfnc(1:idx) + "proc-TwiliteKit-do-make-input-func-nomodel-recalibrated_inputfunc.nii.gz";
                
                if ~isfile(fqfn_twil)
                    continue
                end
                disp(fqfn_twil)
                ic = mlfourd.ImagingContext2(fqfn_twil);
                ic = ic ./ T.inveff(row);
                ic.fqfn = fqfn_recal;
                ic.save()
            end
        end        
        
        function build_all_martin_v1(this)
            import mlkinetics.*

            T = this.table_maframes_nii;
            T = T(T.tracer == "co", :);

            U = this.table_idif_nii();
            U = U(U.tracer == "co", :);

            V = this.table_twilite_nii();
            V = V(V.tracer == "co", :);

            for row = 1:size(T, 1)
                petFqfn = T.bids_fqfn(row);
                U_img = U.img{row};
                V_img = V.img{row};
                this.build_martin_v1(petFqfn, U_img, V_img);
            end
        end
        function build_martin_v1(this, petFqfn, U_img, V_img)
            arguments
                this mlvg.Lee2024
                petFqfn {mustBeFile}
                U_img {mustBeNumeric}
                V_img {mustBeNumeric}
            end

            import mlkinetics.*

            % petDir = ...
            %     fullfile(getenv("SINGULARITY_HOME"), "CCIR_01211", "sourcedata", "sub-108293", "ses-20210421152358", "pet");
            % petFqfn = ...
            %     fullfile(petDir, "sub-108293_ses-20210421152358_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz");

            petMed = mlvg.Ccir1211Mediator.create(petFqfn);

            flirt_out_fqfn = strcat(petMed.fqfp, "-schaeffer.nii.gz");
            if ~isfile(flirt_out_fqfn)
                flirt = mlfsl.Flirt( ...
                    'in', petMed.schaeffer_ic, ...
                    'ref', petMed.imagingReference, ...
                    'out', flirt_out_fqfn, ...
                    'omat', fullfile(petMed.derivPetPath, "T1w_on_" + petMed.imagingReference.fileprefix + ".mat"), ...
                    'bins', 1024, ...
                    'interp', 'nearestneighbour', ...
                    'noclobber', false);
                flirt.applyXfm();
            end

            schBK = BidsKit.create( ...
                bids_tags="ccir1211", bids_fqfn=flirt_out_fqfn);
            rk = RepresentationKit.create( ...
                representation_tags="trivial");
            pk = ParcKit.create( ...
                bids_kit=schBK, representation_kit=rk, parc_tags="schaeffer-schaeffer");

            p = pk.make_parc();
            ic1 = p.reshape_to_parc(petMed.imagingContext);
            %disp(ic1)
            %ic1.view()
            ic1.filepath = strrep(ic1.filepath, "sourcedata", "derivatives");
            if ~isfile(ic1.fqfn)
                ic1.save();
            end

            try
                ifc2 = ic1.imagingFormat;
                % U_img = smoothdata(U_img);
                ifc2.img = 1e-3 * ifc2.img ./ U_img;  % kBq/mL -> Bq/mL
                ifc2.img = mean(ifc2.img(:, 121:end), 2);
                ifc2.fileprefix = ifc2.fileprefix + "-idif_martinv1";
                ifc2.save();
            catch ME
                handwarning(ME)
            end

            try
                ifc3 = ic1.imagingFormat;
                % V_img = smoothdata(V_img);
                Nt = min(size(ifc3.img, 2), length(V_img));
                ifc3.img = 1e-3 * ifc3.img(:, 1:Nt) ./ V_img(1:Nt);  % kBq/mL -> Bq/mL
                ifc3.img = mean(ifc3.img(:, 121:end), 2);
                ifc3.fileprefix = ifc3.fileprefix + "-twilite_martinv1";
                if ~any(isnan(ifc3.img))
                    ifc3.save();
                end
            catch ME
                handwarning(ME)
            end
        end
                
        function build_schaeffer_parc(this, select_tag)
            arguments
                this mlvg.Lee2024
                select_tag {mustBeText}  % "select-brain", "select-gm", ...
            end
            if ~any(startsWith(select_tag, "schaeffer-"))
                select_tag = "schaeffer-" + select_tag;
            end

            import mlkinetics.*

            T = this.table_maframes_timeAppend_nii();
            derivs_path = fullfile(getenv("SINGULARITY_HOME"), "CCIR_01211", "derivatives");
            for row = 1:size(T, 1)

                petdyn_ic = mlfourd.ImagingContext2(T.bids_fqfn(row));

                ic1s = cell(size(select_tag));
                parfor tagi = 1:length(select_tag)
                    parcs = mglob(fullfile(derivs_path, T.sub(row), "*/Parcellations/Schaefer2018_200Parcels_7Networks_order_T1_complete.nii.gz"));
                    assert(~isempty(parcs))
                    bk = BidsKit.create(bids_tags="ccir1211", bids_fqfn=parcs(1));
                    pk = ParcKit.create(bids_kit=bk, parc_tags=select_tag(tagi));
                    p = pk.make_parc();

                    ic1 = p.reshape_to_parc(petdyn_ic);
                    disp(ic1.fqfn)
                    ic1s{tagi} = ic1;
                    %plot(ic1)
                    %save(ic1)
                end
                for tagi = 1:length(select_tag)
                    save(ic1s{tagi});
                end
            end
        end
        function build_schaeffer_4(this)
            mats_path = ...
                fullfile( ...
                getenv("SINGULARITY_HOME"), ...
                "CCIR_01211", "derivatives", "sub-108293", "ses-20210218", "Parcellations");
            mats = mats_path + filesep + ["indices_brain.mat", "indices_gm.mat", "indices_wm.mat", "indices_subcortical.mat"];
            ld = load(mats(1));
            brain = ld.indices_brain;
            ld = load(mats(2));
            gm = ld.indices_gm;
            ld = load(mats(3));
            wm = ld.indices_wm;
            ld = load(mats(4));
            subcortical = ld.indices_subcortical;
            
            sch_path = fullfile(getenv("SINGULARITY_HOME"), "CCIR_01211", "derivatives");
            schaeffers = mglob(fullfile(sch_path, "sub-*", "ses-*", "pet", "sub*_ses*_*-schaeffer-schaeffer.nii.gz"));
            for s = asrow(schaeffers)
                ifc = mlfourd.ImagingFormatContext2(s(1));
                img = nan(4, size(ifc, 2));
                img(1, :) = mean(ifc.img(brain, :), 1);
                img(2, :) = mean(ifc.img(gm, :), 1);
                img(3, :) = mean(ifc.img(wm, :), 1);
                img(4, :) = mean(ifc.img(subcortical, :), 1);
                ifc.img = img;
                ifc.fileprefix = strrep(ifc.fileprefix, "-schaeffer-schaeffer", "-schaeffer-select-4");
                ifc.save();
            end
        end

        function this = Lee2024(varargin)
        end
    end

    methods (Static)
        function ic = append_activity_densities(ic1, ic2)
            %% appends extended activity densities in ic2 onto ic1;
            %  must be consistently decay-corrected or decay-uncorrected;
            %  appends img, timesMid, times, taus;
            %  saves backup of ic1, then overwrites ic1 on the filesystem

            ifc1 = ic1.imagingFormat;
            img1 = ifc1.img;
            j1 = ifc1.json_metadata;
            timesMid1 = j1.timesMid;
            times1 = j1.times;
            taus1 = j1.taus;

            ifc2 = ic2.imagingFormat;
            img2 = ifc2.img;
            j2 = ifc2.json_metadata;
            timesMid2 = j2.timesMid;
            times2 = j2.times;
            taus2 = j2.taus;

            [~,idx_star] = max(timesMid2 > timesMid1(end));
            j1.timesMid = [timesMid1; timesMid2(idx_star:end)];
            j1.times = [times1; times2(idx_star:end)];
            j1.taus = [taus1; taus2(idx_star:end)];
            img1 = [img1, img2(idx_star:end)];
            ifc1.json_metadata = j1;
            ifc1.img = img1;
            ic = mlfourd.ImagingContext2(ifc1);

            % ic1.fileprefix = ic1.fileprefix + "_" + stackstr(use_dashes=true) + "-backup";
            % ic1.save();
            % assert(~strcmp(ifc1.fqfn, ic1.fqfn));  % protect large ic1 which generates ifc1 from handles

            % plot(ic);
            ic.fileprefix = ic.fileprefix + "-embed";
            % ic.save();
        end
        function hct = found_hct(fqfn)
            [~,fp] = myfileparts(fqfn);
            re = regexp(fp, "(?<sub>sub-\d{6})_ses-\d+_\S+", "names");
            switch char(re.sub)
                case 'sub-108293'
                    hct = 46.8;
                case 'sub-108237'
                    hct = 43.9;
                case 'sub-108254'
                    hct = 37.9;
                case 'sub-108250'
                    hct = 42.8;
                case 'sub-108284'
                    hct = 39.7;
                case 'sub-108306'
                    hct = 41.1;
                otherwise
                    error("mlvg:ValueError", stackstr());
            end
        end
        function cifti = qms_to_dtseries(qm_fqfns, new_fqfn, opts)
            arguments
                qm_fqfns {mustBeText}
                new_fqfn {mustBeTextScalar}
                opts.qm_kind {mustBeScalarOrEmpty} = 1
                opts.multiplier single {mustBeScalarOrEmpty} = 1  % for CBF or PS, multiplier ~ 60
                opts.do_write logical = true
            end
            if ~endsWith(new_fqfn, '.dtseries.nii')
                new_fqfn = strcat(new_fqfn, '.dtseries.nii');
            end

            % append cdata from all qm_fqfns
            cifti = mlvg.Lee2024.qm_to_dscalar(qm_fqfns(1), qm_kind=opts.qm_kind, multiplier=opts.multiplier);
            for fni = 2:length(qm_fqfns)
                cifti_ = mlvg.Lee2024.qm_to_dscalar( ...
                    qm_fqfns(fni), qm_kind=opts.qm_kind, multiplier=opts.multiplier);
                cifti.cdata = [cifti.cdata, cifti_.cdata];
            end

            % adjust diminfo{2}
            cifti.diminfo{2} = cifti_diminfo_make_series(size(cifti.cdata, 2), 0, 1, 'SECOND');

            % write
            if opts.do_write
                cifti_write(cifti, convertStringsToChars(new_fqfn));
            end
        end
        function cifti = qm_to_dscalar(qm_fqfn, opts)
            arguments
                qm_fqfn {mustBeFile}
                opts.new_fqfn {mustBeTextScalar} = ""
                opts.qm_kind {mustBeScalarOrEmpty} = 1
                opts.multiplier single {mustBeScalarOrEmpty} = 1  % for CBF or PS, multiplier ~ 60
                opts.do_write logical = false
            end
            if isemptytext(opts.new_fqfn)
                opts.new_fqfn = fullfile( ...
                    myfileparts(qm_fqfn), mybasename(qm_fqfn) + "-kind" + opts.qm_kind + ".dscalar.nii");
            end

            % qm nifti
            ifc = mlfourd.ImagingFormatContext2(qm_fqfn);
            img = ifc.img(:, opts.qm_kind);

            % cifti for supported Schaefer indices
            schaefer = cifti_read( ...
                fullfile( ...
                    getenv("SINGULARITY_HOME"), ...
                    "CCIR_01211", "derivatives", "sub-108293", "ses-20210218", "Parcellations", ...
                    "Schaefer164k.dscalar.nii"));
            schaefer_indices = schaefer.cdata;  % ~ 734132x1 single
            unique_schaefer_indices = unique(schaefer_indices);
            unique_schaefer_indices(unique_schaefer_indices == 0) = [];  % remove zero

            % aufbau new cdata1
            ld = load( ...                
                fullfile( ...
                    getenv("SINGULARITY_HOME"), ...
                    "CCIR_01211", "derivatives", "sub-108293", "ses-20210218", "Parcellations", ...
                    "indices_all.mat"));
            indices_all = ld.indices_all;  % no zero
            cdata1 = zeros(size(schaefer_indices), "single");
            for i_ = asrow(unique_schaefer_indices)
                found = indices_all == i_;
                qm = img(found);  % must be scalar
                cdata1(i_ == schaefer_indices) = qm;
            end
            cdata1 = opts.multiplier * single(cdata1);

            % write
            cifti = schaefer;
            cifti.cdata = cdata1;
            if opts.do_write
                cifti_write(cifti, convertStringsToChars(opts.new_fqfn));
            end
        end
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
        function [nomodelembed_ic,idealembed_ic] = embed_countsFdg_for_TwiliteKit_inputfunc(nomodel_ic, ideal_ic)
            %% manually specified as described in using_countsFdg_for_TwiliteKit-do-make-inputfunc-nomodel.m

            arguments
                nomodel_ic mlfourd.ImagingContext2
                ideal_ic mlfourd.ImagingContext2
            end
            assert(contains(nomodel_ic.fileprefix, "TwiliteKit-do-make-input-func-nomodel_inputfunc"))
            assert(contains(ideal_ic.fileprefix, "TwiliteKit-do-make-input-func-nomodel_inputfunc_dynesty-RadialArtery-ideal"))
            plot(nomodel_ic);
            plot(ideal_ic);

            %% working variables

            halflife = 6586.272;

            %% extend nomodel_ic using CCIRRadMeasurements

            rm_path = fullfile(getenv("HOME"'), "Documents", "CCIRRadMeasurements");

            if contains(nomodel_ic.filepath, "sub-108293")

                rm = mlpet.CCIRRadMeasurements.createFromDate(datetime(2021,4,21));

                %% build tables

                t_drawn = rm.countsFdg.TIMEDRAWN_Hh_mm_ss;
                activity = rm.countsFdg.DECAYCorrSpecificActivity_KBq_mL;  % decay-corrected
                t_elapsed = seconds(t_drawn - datetime(2021,4,21,15,57,9, TimeZone="local")) + 31;
                T = table(t_drawn, t_elapsed, activity);
                T(7:8,:) = [];  % remove empty samples
                T([2,4,6,end],:) = [];  % remove plasma samples
                T.activity = 1e3 * T.activity;  % kBq/mL -> Bq/mL

                %% adjust ifc.img, ifc.json_metadata, ifc.fileprefix
                
                ifc = nomodel_ic.imagingFormat;                

                img = ifc.img;
                %img = img .* 2.^(asrow(ifc.json_metadata.timesMid) / halflife);  % decay-corrected to draw time
                img = [img(1:470), T.activity'];
                ifc.img = img;

                j = ifc.json_metadata;
                j.timesMid = [j.timesMid(1:470); T.t_elapsed];
                j.times = [j.timesMid(1:470); (T.t_elapsed - 1)];
                j.taus = j.taus(1:474);
                ifc.json_metadata = j;

                ifc.fileprefix = ifc.fileprefix + "-embed";
                nomodelembed_ic = mlfourd.ImagingContext2(ifc);

            elseif contains(nomodel_ic.filepath, "sub-108237")

                rm = mlpet.CCIRRadMeasurements.createFromDate(datetime(2022,10,31));

                %% build tables

                t_drawn2 = rm.countsFdgVenous.TIMEDRAWN_Hh_mm_ss;
                activity2 = rm.countsFdgVenous.DECAYCorrSpecificActivity_KBq_mL;  % decay-corrected to draw time
                t_elapsed2 = seconds(t_drawn2 - datetime(2022,10,31,12,8,4)) + 17;
                T2 = table(t_drawn2, t_elapsed2, activity2);

                t_drawn = rm.countsFdg.TIMEDRAWN_Hh_mm_ss;
                activity = rm.countsFdg.DECAYCorrSpecificActivity_KBq_mL;  % decay-corrected
                t_elapsed = seconds(t_drawn - datetime(2022,10,31,12,8,4, TimeZone="local")) + 17;
                T = table(t_drawn, t_elapsed, activity);

                T2.Properties.VariableNames = T.Properties.VariableNames;                
                T.t_drawn.TimeZone = "local";
                T2.t_drawn.TimeZone = "local";

                U = [T; T2];
                U = sortrows(U, "t_elapsed");
                U([2,4,6,8,10,12,14,end],:) = [];  % remove plasma samples
                U.activity = 1e3 * U.activity;  % kBq/mL -> Bq/mL

                %% adjust ifc.img, ifc.json_metadata, ifc.fileprefix

                ifc = nomodel_ic.imagingFormat;

                img = ifc.img;
                %img = img .* 2.^(asrow(ifc.json_metadata.timesMid) / halflife);  % decay-corrected to draw time
                img = [img, asrow(U.activity(2:end))];
                img(393:406) = nan;
                select = ~isnan(img);
                ifc.img = img(select);

                j = ifc.json_metadata;
                j.timesMid = [j.timesMid; U.t_elapsed(2:end)];
                j.times = [j.times; U.t_elapsed(2:end) - 1];
                j.taus = [j.taus; ones(7,1)];
                j.timesMid = j.timesMid(select);
                j.times = j.times(select);
                j.taus = j.taus(select);
                ifc.json_metadata = j;

                ifc.fileprefix = ifc.fileprefix + "-embed";
                nomodelembed_ic = mlfourd.ImagingContext2(ifc);

            elseif contains(nomodel_ic.filepath, "sub-108250")

                rm = mlpet.CCIRRadMeasurements.createFromDate(datetime(2022,12,7));

                %% build tables

                t_drawn2 = rm.countsFdgVenous.TIMEDRAWN_Hh_mm_ss;
                activity2 = rm.countsFdgVenous.DECAYCorrSpecificActivity_KBq_mL;  % decay-corrected to draw time
                t_elapsed2 = seconds(t_drawn2 - datetime(2022,12,7,10,49,9));
                T2 = table(t_drawn2, t_elapsed2, activity2);                
                T2.activity2 = 1e3 * T2.activity2;  % kBq/mL -> Bq/mL
                T2 = T2(1:2:end-1, :);  % remove plasma samples

                %% adjust ifc.img, ifc.json_metadata, ifc.fileprefix

                ifc = nomodel_ic.imagingFormat;

                Tidx = floor(t_elapsed2(1)) - 1;
                img = ifc.img(1:Tidx); 
                %img = img .* 2.^(asrow(ifc.json_metadata.timesMid(1:Tidx)) / halflife);  % decay-corrected to draw time
                img = [img, T2.activity2'];  % extend activity using CCIRRadMeasurements
                ifc.img = img;

                j = ifc.json_metadata;
                j.timesMid = [j.timesMid(1:Tidx); T2.t_elapsed2];
                j.times = [j.times(1:Tidx); T2.t_elapsed2 - 1];
                j.taus = [j.taus(1:Tidx); ones(size(T2.t_elapsed2))];
                ifc.json_metadata = j;

                ifc.fileprefix = ifc.fileprefix + "-embed";
                nomodelembed_ic = mlfourd.ImagingContext2(ifc);

            elseif contains(nomodel_ic.filepath, "sub-108306")

                rm = mlpet.CCIRRadMeasurements.createFromFilename( ...
                    fullfile(rm_path, "CCIRRadMeasurements 20230227.xlsx"));

                %% build tables

                t_drawn2 = rm.countsFdgVenous.TIMEDRAWN_Hh_mm_ss;
                activity2 = rm.countsFdgVenous.DECAYCorrSpecificActivity_KBq_mL;  % decay-corrected to draw time
                t_elapsed2 = seconds(t_drawn2 - datetime(2023,2,27,11,38,53));
                T2 = table(t_drawn2, t_elapsed2, activity2);                
                T2.activity2 = 1e3 * T2.activity2;  % kBq/mL -> Bq/mL
                T2 = T2(1:2:end-1, :);  % remove plasma samples

                %% adjust ifc.img, ifc.json_metadata, ifc.fileprefix

                ifc = nomodel_ic.imagingFormat;

                img = ifc.img; 
                %img = img .* 2.^(asrow(ifc.json_metadata.timesMid) / halflife);  % decay-corrected to draw time
                img = [img, T2.activity2'];  % extend activity using CCIRRadMeasurements
                ifc.img = img;

                j = ifc.json_metadata;
                j.timesMid = [j.timesMid; T2.t_elapsed2];
                j.times = [j.times; T2.t_elapsed2 - 1];
                j.taus = [j.taus; ones(size(T2.t_elapsed2))];
                ifc.json_metadata = j;

                ifc.fileprefix = ifc.fileprefix + "-embed";
                nomodelembed_ic = mlfourd.ImagingContext2(ifc);

            else

                error("mlvg:ValueError", stackstr())

            end
            plot(nomodelembed_ic);
            save(nomodelembed_ic);

            %% extend ideal_ic using nomodelembed_ic

            idealembed_ic = mlvg.Lee2024.append_activity_densities(ideal_ic, nomodelembed_ic); 
            plot(idealembed_ic);
            save(idealembed_ic);
        end

        % Record ID,Event Name,Date of session:,Hematocrit
        % 108237,White Matter Hyperintensities (Arm 4: WMH),10/31/2022,43.9
        % 108250,Measuring Aerobic Glycolysis  (Arm 3: MAG),12/7/2022,42.8
        % 108254,Measuring Aerobic Glycolysis  (Arm 3: MAG),11/16/2022,37.9
        % 108283,White Matter Hyperintensities (Arm 4: WMH),2/15/2023,42.1
        % 108284,White Matter Hyperintensities (Arm 4: WMH),2/20/2023,39.7
        % 108293,Measuring Aerobic Glycolysis  (Arm 3: MAG),4/21/2021,46.8
        % 108306,White Matter Hyperintensities (Arm 4: WMH),2/27/2023,41.1

        function csv2nii(fn)
            csv = readtable(fn);
            sfn = split(fn, "_dynesty");
            nii__ = mlfourd.ImagingFormatContext2(sfn(1) + ".nii.gz");
            j = nii__.json_metadata;
            j1 = j;

            vars = csv.Properties.VariableNames;
            if any(contains(vars, "signal"))
                nii = copy(nii__);
                nii.img = asrow(single(csv.signal)); 
                j1.timesMid = csv.timesMid;
                viable = ~isnan(j.timesMid);
                j1.times = j.times(viable);
                j1.taus = j.taus(viable);
                nii.json_metadata = j1;
                nii.fqfp = myfileprefix(fn);
                nii.save();
            end
            if any(contains(vars, "ideal"))
                nii = copy(nii__);
                nii.img = asrow(single(csv.ideal));
                j1.times = csv.times;
                j1.taus = ones(size(j1.times));
                j1.timesMid = j1.times + j1.taus/2;
                nii.json_metadata = j1;
                nii.fqfp = myfileprefix(fn);
                nii.save();
            end
        end
        function s = remove_tau_tag(s)
            re = regexp(s, "\S+(?<ttag>\-tau\d+)\-\S+", "names");
            s = strrep(s, re.ttag, "");
        end
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
        table_maframes_timeAppend_nii_
        table_twilite_nest_
        table_twilite_nii_

        par_labels__ = ["t_0", "\tau_2", "\tau_3", "\alpha - 1", "1/\beta", "p", "\delta p_2", "\delta p_3", "1/\gamma", "f_2", "f_3", "f_{ss}", "A", "\sigma"]
        par_names__ = "p" + (0:13)
        sub_names__ = ["sub-108293", "sub-108237", "sub-108254", "sub-108250", "sub-108284", "sub-108306"]
        tracers__ = ["co", "oo", "ho", "fdg"]
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
        function fn = fqfn2crv(this, fqfn)
            dt = this.fqfn2datetime(fqfn);
            dstr = string(datetime(dt, Format="yyyyMMdd"));
            trc = this.filename2tracer(mybasename(fqfn));
            switch lower(trc)
                case {'co','oc','ho','oo'}
                    fn = sprintf("o15_dt%s.crv", dstr);
                case 'fdg'
                    fn = sprintf("fdg_dt%s.crv", dstr);
                otherwise
                    error("mlvg:ValueError", "%s: trc=%s not recognized.", stackstr(), trc)
            end
        end
        function T = strrep_bids_fqfn(this, T, opts)
            %% updates this.table_maframes_nii_, replacing arbitrary strings with new strings

            arguments
                this mlvg.Lee2024
                T table
                opts.s1 {mustBeTextScalar} = fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity")
                opts.s2 {mustBeTextScalar} = fullfile(getenv("SINGULARITY_HOME"))
            end
            
            T.bids_fqfn = strrep(T.bids_fqfn, opts.s1, opts.s2);
        end
    end

    methods (Static) %, Access = private)
        function T = add_csv(T, do_rescale)
            arguments
                T table
                do_rescale logical = false  % use T.inveff to rescale img
            end

            timesMid = cell(size(T, 1), 1);
            img = cell(size(T, 1), 1);
            for row = 1:size(T, 1)
                csv = readtable(T.bids_fqfn(row));
                try
                    timesMid{row} = double(asrow(csv.times));
                catch ME
                    handwarning(ME)
                    timesMid{row} = double(asrow(csv.timesMid));
                end
                try
                    img__ = double(asrow(csv.ideal));
                catch ME
                    handwarning(ME)
                    img__ = double(asrow(csv.img));
                end
                if do_rescale
                    inveff__ = T.inveff(row);
                    img{row} = img__ .* inveff__ / 1e3; % kBq/mL    
                else
                    img{row} = img__ / 1e3; % kBq/mL    
                end
            end
            T = addvars(T, timesMid, img, NewVariableNames={'timesMid', 'img'});
        end
        function T = add_imaging(T, do_rescale)
            arguments
                T table
                do_rescale logical = false  % use T.inveff to rescale img
            end

            imaging = cell(size(T, 1), 1);
            timesMid = cell(size(T, 1), 1);
            img = cell(size(T, 1), 1);
            for row = 1:size(T, 1)
                ic = mlfourd.ImagingContext2(T.bids_fqfn(row));
                ic.selectImagingTool();
                imaging{row} = ic;
                timesMid__ = asrow(ic.json_metadata.timesMid);
                img__ = asrow(double(ic));
                viable = ~isnan(double(timesMid__)); % & (img_ > 1); % noise ~ 1 kBq/mL
                timesMid{row} = timesMid__(viable);
                if do_rescale
                    inveff_ = T.inveff(row);
                    img{row}  = img__(viable) .* inveff_ / 1e3; % kBq/mL
                else
                    img{row}  = img__(viable) / 1e3; % kBq/mL
                end
            end
            T = addvars(T, imaging, timesMid, img, NewVariableNames={'imaging', 'timesMid', 'img'});
        end
        function T = add_peak_position(T)
            arguments
                T table
            end

            pp = nan(size(T, 1), 1);
            for row = 1:size(T, 1)
                [~, pp(row)] = max(T.img{row});
                if pp(row) == length(T.img{row})
                    pp(row) = 1;
                end
            end
            T = addvars(T, pp, NewVariableNames={'peak_position'});
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
        function dt = fqfn2datetime(fqfn)
            arguments
                fqfn {mustBeTextScalar}
            end

            parts = split(fqfn, filesep);
            selected = contains(parts, "ses-");
            [~,sidx] = max(selected);
            selected(sidx+1:end) = false;
            ses = parts(selected);
            dt = datetime(extractBetween(ses, 5, 18), InputFormat="yyyyMMddHHmmss");
        end
        function N = filename2Nkind(fn)
            [~,fp] = myfileparts(fn);
            if contains(fp, "Raichle") || contains(fp, "Mintun")
                N = 6; return
            end
            if contains(fp, "Mintun")
                N = 7; return
            end
            if contains(fp, "Ichise")
                N = 9; return
            end
            error("mlvg:ValueError", stackstr())
        end
        function trc = filename2tracer(fn)
            arguments
                fn {mustBeTextScalar}
            end

            try
                re = regexp(fn, "\S+_trc-(?<trc>[a-zA-Z0-9]+)_\S+", "names");
                trc = re.trc;
            catch ME
                handwarning(ME)
                trc = "unknown";
            end
        end
        function m = infer_multiplier(fqfn, kind)
            [~,fp] = myfileparts(fqfn);
            if contains(fp, "Raichle")
                switch kind
                    case {1,3}
                        m = 60; return
                    otherwise
                        m = 1; return
                end
            end
            m = 1;
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
        function lims = xlim(trc)
            switch char(trc)
                case {'co', 'oc'}
                    lims = [-80, 120];
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
        function lims = ylim(trc, is_ideal)
            arguments
                trc {mustBeTextScalar}
                is_ideal logical = false
            end

            if is_ideal
                switch char(trc)
                    case {'co', 'oc'}
                        lims = [0, 1000];
                    case 'oo'
                        lims = [0, 1000];
                    case 'ho'
                        lims = [0, 1000];
                    case 'fdg'
                        lims = [0, 380];

                    otherwise
                        error("mlvg:RuntimeError", ...
                            "%s: %s", stackstr(), trc)
                end
            else
                switch char(trc)
                    case {'co', 'oc'}
                        lims = [0, 500];
                    case 'oo'
                        lims = [0, 500];
                    case 'ho'
                        lims = [0, 500];
                    case 'fdg'
                        lims = [0, 190];

                    otherwise
                        error("mlvg:RuntimeError", ...
                            "%s: %s", stackstr(), trc)
                end
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
