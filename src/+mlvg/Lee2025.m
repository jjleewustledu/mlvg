classdef Lee2025 < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 04-Jul-2025 01:08:28 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.
   

    properties (Constant)
        % PARC_SCHAEF_TAG = "-ParcSchaeffer-reshape-to-schaeffer-schaeffer"
        PARC_SCHAEF_TAG = "-ParcSchaeffer-invariant-schaeffer-schaeffer"
        % PARC_SCHAEF_TAG = "-ParcSchaeffer-highsnr-schaeffer-schaeffer"
    end

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
                opts.steps {mustBeNumericOrLogical} = 5
            end

            % exclusions
            if 5 == opts.steps
                fn = extractBefore(mybasename(opts.nii), "_proc") + "_proc-MipIdif_idif.nii.gz";
                target = fullfile(fileparts(opts.nii), fn);
                target = strrep(target, "sourcedata", "derivatives");
                if isfile(target)
                    fprintf("%s: skipping existing %s\n", stackstr(), target);
                    return
                end
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
        function durations = build_schaeffer_parc(fqfns, opts)
            %% e.g.,
            %  sub-108306_ses-20230227113853_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz ->
            %  sub-108306_ses-20230227113853_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames-ParcSchaeffer-reshape-to-schaeffer-schaeffer.nii.gz

            arguments
                fqfns {mustBeText}
                opts.out_dir {mustBeFolder} = fullfile(getenv("SINGULARITY_HOME"), "CCIR_01211")
                opts.do_plot logical = false
            end
            if ~contains(fqfns(1), opts.out_dir)
                fqfns = fullfile(opts.out_dir, fqfns);
            end

            import mlkinetics.*

            durations = nan(1, length(fqfns));

            for fidx = 1:length(fqfns)
                try
                    tic

                    fqfn = fqfns(fidx);
                    petMed = mlvg.Ccir1211Mediator.create(fqfn);

                    trc = lower(petMed.tracer);
                    foundT1w = mglob(fullfile(petMed.derivPetPath, sprintf("T1w_on_*trc-%s*.nii.gz", trc)));
                    if isempty(foundT1w)
                        continue
                    end
                    foundRef = mglob(fullfile(petMed.sourcePetPath, extractAfter(mybasename(foundT1w, withext=true), "T1w_on_")));
                    if isempty(foundRef)
                        continue
                    end
                    imagingReference = mlfourd.ImagingContext2(foundRef);
                    schaef_flirted_fqfn = strcat(petMed.fqfp, "-schaeffer.nii.gz");
                    schaef_flirted_fqfn = strrep(schaef_flirted_fqfn, "sourcedata", "derivatives");
                    target_fqfn = strrep(schaef_flirted_fqfn, "-schaeffer", mlvg.Lee2025.PARC_SCHAEF_TAG);
                    if isfile(target_fqfn)
                        continue
                    end
                    omat = fullfile(petMed.derivPetPath, mybasename(foundT1w) + ".mat");
                    if ~isfile(omat)
                        omat = mglob(fullfile(petMed.derivSubPath, "ses-*", "pet", sprintf("T1w_on_%s.mat", trc)));
                    end
                    assert(isfile(omat))
                    flirt = mlfsl.Flirt( ...
                        'in', petMed.schaeffer_ic, ...
                        'ref', imagingReference, ...
                        'out', schaef_flirted_fqfn, ...
                        'omat', omat, ...
                        'bins', 1024, ...
                        'interp', 'nearestneighbour', ...
                        'noclobber', true);
                    flirt.applyXfm();

                    bk = BidsKit.create(bids_tags="ccir1211", bids_fqfn=schaef_flirted_fqfn);
                    pk = ParcKit.create(bids_kit=bk, parc_tags="schaeffer-schaeffer");
                    p = pk.make_parc();

                    ic1 = p.reshape_to_parc_fast(fqfn);  % petMed.imagingContext                    
                    ic1.save();
                    
                    if opts.do_plot
                        disp(ic1.fqfn)
                        plot(ic1')
                    end

                    durations(fidx) = toc;
                catch ME
                    handwarning(ME)
                end
            end
        end

        function ic = build_schaeffer_finite(nmaf, opts)
            %% use with build_schaeffer_parc

            arguments
                nmaf {mustBeText}
                opts.noclobber logical = true
            end

            ic = [];

            try
                % mark tracers with multiple delay sets
                if contains(nmaf, "trc-oo") || contains(nmaf, "trc-fdg")
                    nmaf = regexprep(nmaf, "-delay\d+", "-delay*");
                end

                % check noclobber
                [pth,fp] = myfileparts(nmaf);
                pth_derivs = strrep(pth, "sourcedata", "derivatives");
                if contains(fp, "-delay")
                    fp_final = extractBefore(fp, "-delay") + mlvg.Lee2025.PARC_SCHAEF_TAG + "-finite";
                else
                    fp_final = extractBefore(fp, "proc-") + mlvg.Lee2025.PARC_SCHAEF_TAG + "-finite";
                end
                fqfn_final = fullfile(pth_derivs, fp_final + ".nii.gz");
                if opts.noclobber && isfile(fqfn_final)
                    ic = mlfourd.ImagingContext2(fqfn_final);
                    return
                end

                if isfile(nmaf)
                    ic = base_case(nmaf);
                    ic.fqfn = fqfn_final;
                    ic.save();
                    return
                end
                if contains(nmaf, "*")
                    nmafs = mglob(nmaf);
                    nmafs = nmafs(~contains(nmafs, "_avgt") & ~contains(nmafs, "_mipt"));
                    nmafs = natsort(nmafs);
                    ics = [];
                    for n = nmafs
                        ics = [ics, base_case(n)]; %#ok<AGROW>
                    end
                    ic = ics(1);
                    for iidx = 2:length(nmafs)
                        ic = ic.timeAppend(ics(iidx), concat_json_fields=true);
                    end
                    ic.fqfn = fqfn_final;
                    ic.save();
                    return
                end
            catch ME
                handwarning(ME)
            end

            function ic = base_case(nmaf)
                ic = [];
                try
                    [~,fp_] = myfileparts(nmaf);
                    if contains(fp_, "-BrainMoCo2")
                        fp_ = extractBefore(fp_, "-BrainMoCo2") + "-BrainMoCo2-createNiftiMovingAvgFrames" + mlvg.Lee2025.PARC_SCHAEF_TAG;  % retain "-delay\d+"
                    elseif contains(fp_, "-console")
                        fp_ = extractBefore(fp_, "-console") + "-consoleDynamic" + mlvg.Lee2025.PARC_SCHAEF_TAG;
                    else
                        error("mlvg:ValueError", stackstr())
                    end
                    fqfn_ = fullfile(pth_derivs, fp_ + ".nii.gz");

                    ic = mlfourd.ImagingContext2(fqfn_);
                    ic = mlpipeline.ImagingMediator.ensureFiniteImagingContext2(ic);  % reformat json, remove null|nan
                catch ME
                    handwarning(ME);
                end
            end
        end

        function [ic,ic1,ff] = build_nonzero_2d(obj, opts)
            %% scrubs zeros, but allows earliest frames to be zero.
            %  Supply externally determineed faulty_frames as 1-form.
            %  Supply paired obj with identical timings to have identical frame censoring.
            %  For exploration of faults, set do_only_count_faulty = true.
            %  uncensored = 1:20 => no censoring for frames 1:20.
            %  parc_tol defines parcel activity that is nearly zero.
            %  tol is the fraction of zero-parcels in frame to tolerate

            arguments
                obj {mustBeNonempty}
                opts.faulty_frames = []
                opts.paired_obj = []
                opts.do_backup logical = true
                opts.uncensored {mustBeNumeric} = 1:20
                opts.parc_tol {mustBeScalarOrEmpty} = 100  % Bq/cc
                opts.tol {mustBeScalarOrEmpty} = 0.05
                opts.do_only_count_faulty logical = false
            end

            % find zeros
            try
                ic = mlfourd.ImagingContext2(obj);
                img = ic.imagingFormat.img;
                assert(ismatrix(img))
                if ~isempty(opts.faulty_frames)
                    faulty_frames = opts.faulty_frames;
                else
                    Ngo = size(img, 1);  % num. greyordinates
                    iszero2d = abs(img) < opts.parc_tol;
                    iszero2d(:, opts.uncensored) = false;  % be permissive with earliest frames
                    if ~any(iszero2d(2:end))
                        ic1 = []; ff = 0;
                        return
                    end
                    faulty_frames = sum(iszero2d, 1) > opts.tol*Ngo;  % 1d row
                end
                ff = sum(faulty_frames);
                if opts.do_only_count_faulty
                    ic1 = [];
                    return
                end
            catch ME
                ic = []; ic1 = []; ff = nan;
                if isfile(obj)
                    id = obj;
                else
                    id = class(obj);
                end
                fprintf("%s: %s\n", id, ME.message)
                return
            end            

            % backup as requested
            if opts.do_backup
                ic_bak = copy(ic);
                ic_bak.fileprefix = ic.fileprefix + "_bak";
                idx = 1;
                while isfile(ic_bak.fqfn)
                    idx = idx + 1;
                    ic_bak.fileprefix = ic.fileprefix + "_bk" + idx;
                end
                ic_bak.save();
            end
            
            % build nonzero
            img = img(:, ~faulty_frames);
            ifc = ic.imagingFormat;
            ifc.img = img;
            ic = mlfourd.ImagingContext2(ifc);
            j = ic.json_metadata;
            try
                j.starts = j.starts(ascol(~faulty_frames));
            catch ME
                handwarning(ME)
            end
            try
                j.taus = j.taus(ascol(~faulty_frames));
            catch ME
                handwarning(ME)
            end
            try
                j.times = j.times(ascol(~faulty_frames));
            catch ME
                handwarning(ME)
            end
            try
                j.timesMid = j.timesMid(ascol(~faulty_frames));
            catch ME
                handwarning(ME)
            end
            ic.json_metadata = j;

            % save built
            ic.save();

            % build nonzero for paired obj
            if ~isempty(opts.paired_obj)
                try
                    ic1 = mlvg.Lee2025.build_nonzero_2d( ...
                        opts.paired_obj, faulty_frames=faulty_frames, do_backup=opts.do_backup, tol=opts.tol);
                catch ME
                    handwarning(ME)
                end
            else
                ic1 =[];
            end
        end

        function ic = build_mip_idif_finite(nmaf, opts)
            arguments
                nmaf {mustBeText}
                opts.noclobber logical = true
            end

            ic = [];

            try
                % mark tracers with multiple delay sets
                if contains(nmaf, "trc-oo") || contains(nmaf, "trc-fdg")
                    nmaf = regexprep(nmaf, "-delay\d+", "-delay*");
                end

                % check noclobber
                [pth,fp] = myfileparts(nmaf);
                pth_derivs = strrep(pth, "sourcedata", "derivatives");
                if contains(fp, "-delay")
                    fp_final = extractBefore(fp, "-delay") + "-MipIdif-finite_idif";
                else
                    fp_final = extractBefore(fp, "proc-") + "proc-MipIdif-finite_idif";
                end
                fqfn_final = fullfile(pth_derivs, fp_final + ".nii.gz");
                if opts.noclobber && isfile(fqfn_final)
                    ic = mlfourd.ImagingContext2(fqfn_final);
                    return
                end

                if isfile(nmaf)
                    ic = base_case(nmaf);
                    ic.fqfn = fqfn_final;
                    ic.save();
                    return
                end
                if contains(nmaf, "*")
                    nmafs = mglob(nmaf);
                    nmafs = nmafs(~contains(nmafs, "_avgt") & ~contains(nmafs, "_mipt"));
                    nmafs = natsort(nmafs);
                    ics = [];
                    for n = nmafs
                        ics = [ics, base_case(n)]; %#ok<AGROW>
                    end
                    ic = ics(1);
                    for iidx = 2:length(nmafs)
                        ic = ic.timeAppend(ics(iidx), concat_json_fields=true);
                    end
                    ic.fqfn = fqfn_final;
                    ic.save();
                    return
                end
            catch ME
                handwarning(ME)
            end

            function ic = base_case(nmaf)
                ic = [];
                try
                    % load NIfTI moving avg frames
                    assert(isfile(nmaf), "missing nmaf")
                    nmaf_ifc = mlfourd.ImagingFormatContext2(nmaf);

                    % load centerline_on_pet
                    cl_fqfn = fullfile(pth_derivs, "centerline_on_pet.nii.gz");
                    assert(isfile(cl_fqfn), "missing centerline")
                    try
                        cl_ifc = mlfourd.ImagingFormatContext2(cl_fqfn);
                        cl_img = cl_ifc.img;
                    catch ME
                        fprintf("%s: trouble opening fsl edits to nifti: %s\n", stackstr(), ME.message);
                        cl_img = niftiread(cl_fqfn);
                    end

                    % sample centerline from nmaf
                    sz = size(nmaf_ifc.img);
                    nmaf_mat = reshape(nmaf_ifc.img, [sz(1) * sz(2) * sz(3), sz(4)]);
                    cl_vec = reshape(cl_img, [sz(1) * sz(2) * sz(3), 1]);
                    nmaf_mat = nmaf_mat(logical(cl_vec), :);
                    nmaf_form = mean(nmaf_mat, 1);
                    nmaf_ifc.img = nmaf_form;  % reuse malloc
                    [~,fp_] = myfileparts(nmaf);
                    nmaf_ifc.filepath = pth_derivs;
                    if contains(fp_, "-BrainMoCo2")
                        nmaf_ifc.fileprefix = extractBefore(fp_, "-BrainMoCo2") + "-MipIdif_idif";  % retain "-delay\d+"
                    elseif contains(fp_, "-console")
                        nmaf_ifc.fileprefix = extractBefore(fp_, "-console") + "-MipIdif_idif";
                    else
                        error("mlvg:ValueError", stackstr())
                    end
                    ic = mlfourd.ImagingContext2(nmaf_ifc);
                    ic = mlpipeline.ImagingMediator.ensureFiniteImagingContext2(ic);  % reformat json, remove null|nan
                catch ME
                    handwarning(ME);
                end
            end
        end

        function build_all_martin_v1_idif(mat_file)
            arguments
                mat_file {mustBeFile} = fullfile( ...
                    getenv("SINGULARITY_HOME"), "CCIR_01211", "mlvg_Lee2025Par_globbing_co.mat")
            end

            import mlkinetics.*

            ld = load(mat_file);
            globbed_co = ld.globbed;

            for g = globbed_co
                try
                    petFqfn = g;
                    [pth,fp] = myfileparts(petFqfn);
                    pth = strrep(pth, "sourcedata", "derivatives");
                    fp = extractBefore(fp, "-delay0-BrainMoCo2-createNiftiMovingAvgFrames");

                    petSchaeferFqfn = fullfile(pth, fp + mlvg.Lee2025.PARC_SCHAEF_TAG + "-finite.nii.gz");
                    mipidifFqfn = fullfile(pth, fp + "-MipIdif-finite_idif.nii.gz");
                    mlvg.Lee2025.build_martin_v1(petSchaeferFqfn, mipidifFqfn);
                catch ME
                    handwarning(ME)
                end
            end
        end

        function ic = build_martin_v1(petSchaeferFqfn, mipidifFqfn)
            arguments
                petSchaeferFqfn {mustBeFile}
                mipidifFqfn {mustBeFile}
            end

            import mlkinetics.*

            ic_parc = mlfourd.ImagingContext2(petSchaeferFqfn);
            ic_idif = mlfourd.ImagingContext2(mipidifFqfn);

            try
                ifc_parc = ic_parc.imagingFormat;
                timesMid = ifc_parc.json_metadata.timesMid;
                tidx = find(timesMid >= 120);  % use timesMid > 120 sec for equilibration
                tidx = tidx(1);
                parc_img = mean(ifc_parc.img(:, tidx:end), 2, "omitnan");                

                ifc_idif = ic_idif.imagingFormat;
                timesMid = ifc_idif.json_metadata.timesMid;
                tidx = find(timesMid >= 120);  % use timesMid > 120 sec for equilibration
                tidx = tidx(1);
                idif_img = mean(ifc_idif.img(tidx:end));

                % build v1

                ifc_parc.img = parc_img ./ idif_img;
                ifc_parc.fileprefix = ifc_parc.fileprefix + "-idif_martinv1";
                ifc_parc.save();

                % return ic
                ic = mlfourd.ImagingContext2(ifc_parc);
            catch ME
                handwarning(ME)
            end
        end

        function ensure_finite(nmaf_fqfns, opts)
            arguments
                nmaf_fqfns {mustBeText}
                opts.target {mustBeTextScalar} = "schaefer"  % "mipidif"
            end

            for fqfn = nmaf_fqfns
                try
                    switch convertStringsToChars(lower(opts.target))
                        case 'mipidif'
                            [pth,fp] = myfileparts(fqfn);
                            pth = strrep(pth, "sourcedata", "derivatives");
                            fp = extractBefore(fp, "-delay") + "-MipIdif_idif";
                            fqfn1 = fullfile(pth, fp + ".nii.gz");
                            if ~isfile(fqfn1)
                                continue  % no source file
                            end
                            if isfile(strrep(fqfn1, "-MipIdif", "-MidIdif-finite"))
                                continue  % target file already exists
                            end
                        case 'schaefer'
                            [pth,fp] = myfileparts(fqfn);
                            pth = strrep(pth, "sourcedata", "derivatives");
                            fp = extractBefore(fp, "-delay") + "-BrainMoCo2-createNiftiMovingAvgFrames" + mlvg.Lee2025.PARC_SCHAEF_TAG;
                            fqfn1 = fullfile(pth, fp + ".nii.gz");
                            if ~isfile(fqfn1)
                                continue  % no source file
                            end
                            if isfile(strrep(fqfn1, "-schaeffer-schaeffer", "-schaeffer-schaeffer-finite"))
                                continue  % target file already exists
                            end
                        otherwise
                            continue
                    end
                    ic = mlfourd.ImagingContext2(fqfn1);
                    ic = mlpipeline.ImagingMediator.ensureFiniteImagingContext2(ic);
                    if ~isfile(ic.fqfn)
                        ic.save();
                    end
                catch ME
                    handwarning(ME)
                end
            end
        end

        function tmp1 = ensure_martinv1(martinv1)
            arguments
                martinv1 = []
            end

            if isempty(martinv1)
                martinv1 = mglob("derivatives/sub-*/ses-*/pet/sub-*_ses-*_trc-co_*" + ...
                    mlvg.Lee2025.PARC_SCHAEF_TAG + ...
                    "-finite-idif_martinv1.nii.gz");
            end
            tmp1 = [];
            for m = martinv1
                ic = mlfourd.ImagingContext2(m);
                mean_v1 = dipmean(ic);
                if mean_v1 > 0.5
                    fprintf("%s:  mean(v1) -> %g\n", m, mean_v1);
                    plot(ic);
                    tmp1 = [tmp1, m]; 
                end
            end
        end

        function ensure_nonzero_2d(parcschaef_co, mipidif_co, opts)
            %% e.g., for co

            arguments
                parcschaef_co {mustBeText}
                mipidif_co {mustBeText}
                opts.parc_tol {mustBeScalarOrEmpty} = 1e3
            end            
            
            badframes_co = [];
            for p = parcschaef_co
                try
                    [~,~,ff] = mlvg.Lee2025.build_nonzero_2d(p, do_only_count_faulty=true, parc_tol=opts.parc_tol);
                    if ~isnan(ff) && ff > 0; badframes_co = [badframes_co, p]; end
                catch ME
                    fprintf("%s: %s\n", p, ME.message);
                end
            end
            disp(badframes_co)

            for idx = 1:length(badframes_co)
                mlvg.Lee2025.build_nonzero_2d(badframes_co(idx), parc_tol=opts.parc_tol, paired_obj=mipidif_co); 
            end

            badframes_co = [];
            for p = parcschaef_co
                try
                    [~,~,ff] = mlvg.Lee2025.build_nonzero_2d(p, do_only_count_faulty=true, parc_tol=opts.parc_tol);
                    if ~isnan(ff) && ff > 0; badframes_co = [badframes_co, p]; end
                catch ME
                    fprintf("%s: %s\n", p, ME.message);
                end
            end
            disp(badframes_co)
        end

        function fqfn = find_fdg_mipt(other_nii)
            %% e.g., sub-108334_ses-20241216113854_trc-fdg_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_mipt.nii.gz

            arguments
                other_nii {mustBeFile}
            end

            other_nii = strrep(other_nii, "sourcedata", "derivatives");
            derivs_subpth = extractBefore(fileparts(other_nii), filesep + "ses-");
            src_subpth = strrep(derivs_subpth, "derivatives", "sourcedata");
            globbed = mglob(fullfile( ...
                derivs_subpth, "ses-*", "pet", "sub-*_ses-*_trc-fdg_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_mipt.nii.gz"));
            if isscalar(globbed) && ~isemptytext(globbed)
                fqfn = globbed;
                return
            end

            globbed = mglob(fullfile( ...
                src_subpth, "ses-*", "pet", "sub-*_ses-*_trc-fdg_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz"));
            assert(~isempty(globbed))
            assert(isscalar(globbed))
            ic = mlfourd.ImagingContext2(globbed);
            ic = max(ic, [], 4);  % mipt
            ic.fileprefix = strrep(ic.fileprefix, "_max4", "_mipt");
            ic.relocateToDerivativesFolder();
            ic.save();
            fqfn = ic.fqfn;
        end

        function fqfn = find_ho_avgt(other_nii, opts)
            %% e.g., sub-108034_ses-20230717140707_trc-ho_proc-delay0-BrainMoCo2-createNiftiStatic.nii.gz

            arguments
                other_nii {mustBeFile}
                opts.selection {mustBeInteger} = 1
            end

            subpth = extractBefore(fileparts(other_nii), filesep + "ses-");
            globbed = mglob(fullfile(subpth, "ses-*", "pet", "sub-*_ses-*_trc-ho_proc-delay0-BrainMoCo2-createNiftiStatic.nii.gz"));
            assert(~isempty(globbed))
            assert(isscalar(globbed(opts.selection)))
            fqfn = globbed(opts.selection);
        end

        function fqfn = find_t1w(other_nii)
            %% e.g., sub-108333_ses-20241122094951_T1w_MPR_vNav_4e_RMS_orient-std.nii.gz

            arguments
                other_nii {mustBeFile}
            end

            subpth = extractBefore(fileparts(other_nii), filesep + "ses-");
            subpth = strrep(subpth, "sourcedata", "derivatives");
            globbed = mglob(fullfile(subpth, "ses-*", "anat", "sub-*_ses-*_T1w_MPR_vNav_4e_RMS_orient-std.nii.gz"));
            assert(~isempty(globbed))
            assert(isscalar(globbed))
            fqfn = globbed;
        end

        function fqfn = find_t1w_on_fdg(other_nii)
            %% e.g., T1w_on_sub-108333_ses-20241202_trc-fdg_proc-consoleDynamic_avgt.nii.gz

            arguments
                other_nii {mustBeFile}
            end

            subpth = extractBefore(fileparts(other_nii), filesep + "ses-");
            derivs_subpth = strrep(subpth, "sourcedata", "derivatives");
            globbed = mglob(fullfile(derivs_subpth, "ses-*", "pet", "T1w_on_*_trc-fdg_*.nii.gz"));
            assert(~isempty(globbed))
            fqfn = globbed(end);  % console and e7 fdg (ses-\d{14}) may co-exist
        end

        function fqfn = find_t1w_on_ho(other_nii, opts)
            %% e.g., T1w_on_sub-108034_ses-20230717140707_trc-ho_proc-delay0-BrainMoCo2-createNiftiStatic.nii.gz

            arguments
                other_nii {mustBeFile}
                opts.selection {mustBeInteger} = 1
            end

            subpth = extractBefore(fileparts(other_nii), filesep + "ses-");
            subpth = strrep(subpth, "sourcedata", "derivatives");
            globbed = mglob(fullfile(subpth, "ses-*", "pet", "T1w_on_*_trc-ho_proc-delay0-BrainMoCo2-createNiftiStatic.nii.gz"));
            assert(~isempty(globbed))
            assert(isscalar(globbed(opts.selection)))
            fqfn = globbed(opts.selection);
        end

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

        function flirt_t1w(fqfn, opts)
            %% intended for failures of MipIdif;
            %  try again using ho_static -> oo_delay30_avgt, then applyXfm(T1w);
            %  e.g., sub-108034_ses-20230717133405_trc-oo_proc-delay30-BrainMoCo2-createNiftiStatic.nii.gz

            arguments
                fqfn {mustBeFile}
                opts.specialize_for_tracer logical = true
            end

            if contains(fqfn, "sub-108259_ses-20230731133714")
                % OO administered but scanning performed with CO protocol
                opts.specialize_for_tracer = false;
            end
            if opts.specialize_for_tracer
                if contains(fqfn, "trc-co")
                    mlvg.Lee2025.flirt_t1w_on_co(fqfn);
                    return
                end
                if contains(fqfn, "trc-oo")
                    mlvg.Lee2025.flirt_t1w_on_oo(fqfn);
                    return
                end
            end

            import mlvg.Lee2025.mat

            fqfn = strrep(fqfn, "createNiftiMovingAvgFrames", "createNiftiStatic");

            % flirt fdg_mipt -> co_mipt
            t1w_fqfn = mlvg.Lee2025.find_t1w(fqfn);
            [pth,fp] = myfileparts(fqfn);
            pth_derivs = strrep(pth, "sourcedata", "derivatives");
            t1w_on_tracer = fullfile(pth_derivs, "T1w_on_" + fp + ".nii.gz");
            flirt = mlfsl.Flirt( ...
                'in', t1w_fqfn, ...
                'ref', fqfn, ...
                'out', t1w_on_tracer, ...
                'omat', mat(t1w_on_tracer), ...
                'bins', 1024, ...
                'cost', 'normmi', ...
                'dof', 6, ...
                'interp', 'spline', ...
                'noclobber', false);
            ensuredir(fileparts(t1w_on_tracer))
            if ~isfile(t1w_on_tracer)
                flirt.flirt();
                assert(isfile(t1w_on_tracer))
            end
        end

        function flirt_t1w_on_co(co_fqfn, opts)
            %% intended for failures of MipIdif;
            %  try again using fdg_mipt -> co_mipt, then applyXfm(T1w)

            arguments
                co_fqfn {mustBeFile}
                opts.noclobber logical = true
            end

            co_fqfn = strrep(co_fqfn, "createNiftiMovingAvgFrames", "createNiftiStatic");
            co_fqfn = strrep(co_fqfn, "consoleDynamic", "consoleStatic");

            import mlvg.Lee2025.find_fdg_mipt
            import mlvg.Lee2025.find_t1w
            import mlvg.Lee2025.find_t1w_on_fdg
            import mlvg.Lee2025.mat

            fdg_mipt = find_fdg_mipt(co_fqfn);
            fdg_mipt_on_co = myfileprefix(fdg_mipt) + "_on_co.nii.gz";

            % flirt fdg_mipt -> co_mipt
            flirt = mlfsl.Flirt( ...
                'in', fdg_mipt, ...
                'ref', co_fqfn, ...
                'out', fdg_mipt_on_co, ...
                'omat', mat(fdg_mipt_on_co), ...
                'bins', 1024, ...
                'cost', 'normmi', ...
                'dof', 6, ...
                'interp', 'spline', ...
                'noclobber', false);
            if ~opts.noclobber || ~isfile(fdg_mipt_on_co)
                flirt.flirt();
            end
            assert(isfile(mat(fdg_mipt_on_co)))

            % concatXfm & applyXfm
            t1w2fdg_xfm = mat(find_t1w_on_fdg(fdg_mipt));
            src_pth = fileparts(co_fqfn);
            derivs_pth = strrep(src_pth, "sourcedata", "derivatives");
            t1w_on_co = fullfile(derivs_pth, "T1w_on_" + mybasename(co_fqfn, withext=true));
            if isfile(t1w_on_co)
                return
            end
            flirt.concatXfm(AtoB=t1w2fdg_xfm);
            flirt.in = find_t1w(co_fqfn);
            flirt.ref = co_fqfn;
            flirt.out = t1w_on_co;
            ensuredir(fileparts(t1w_on_co))
            flirt.applyXfm();
            assert(isfile(t1w_on_co))
        end

        function flirt_t1w_on_oo(oo_fqfn, opts)
            %% intended for failures of MipIdif;
            %  try again using ho_static -> oo_delay30_avgt, then applyXfm(T1w);
            %  e.g., sub-108034_ses-20230717133405_trc-oo_proc-delay30-BrainMoCo2-createNiftiStatic.nii.gz

            arguments
                oo_fqfn {mustBeFile}
                opts.noclobber logical = false
            end

            import mlvg.Lee2025.find_ho_avgt
            import mlvg.Lee2025.find_t1w
            import mlvg.Lee2025.find_t1w_on_ho
            import mlvg.Lee2025.mat

            oo_fqfn = strrep(oo_fqfn, "delay0", "delay30");
            oo_fqfn = strrep(oo_fqfn, "createNiftiMovingAvgFrames", "createNiftiStatic");

            ho_avgt = mlvg.Lee2025.find_ho_avgt(oo_fqfn);
            ho_avgt_on_oo = myfileprefix(ho_avgt) + "_on_oo.nii.gz";
            ho_avgt_on_oo = strrep(ho_avgt_on_oo, "sourcedata", "derivatives");

            % flirt ho_avgt -> oo
            flirt = mlfsl.Flirt( ...
                'in', ho_avgt, ...
                'ref', oo_fqfn, ...
                'out', ho_avgt_on_oo, ...
                'omat', mat(ho_avgt_on_oo), ...
                'bins', 1024, ...
                'cost', 'normmi', ...
                'dof', 6, ...
                'interp', 'spline', ...
                'noclobber', false);
            if ~opts.noclobber || ~isfile(ho_avgt_on_oo)
                flirt.flirt();
                assert(isfile(ho_avgt_on_oo))
            end

            % concatXfm & applyXfm
            t1w2ho_xfm = mat(find_t1w_on_ho(ho_avgt));
            pth = fileparts(oo_fqfn);
            pth = strrep(pth, "sourcedata", "derivatives");
            t1w_on_oo = fullfile(pth, "T1w_on_" + mybasename(oo_fqfn, withext=true));
            flirt.concatXfm(AtoB=t1w2ho_xfm);
            flirt.in = find_t1w(oo_fqfn);
            flirt.ref = oo_fqfn;
            flirt.out = t1w_on_oo;
            ensuredir(fileparts(t1w_on_oo))
            flirt.applyXfm();
            assert(isfile(t1w_on_oo))
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

        function fqfn = mat(obj)
            if isfile(obj)
                fqfn = strrep(obj, ".nii.gz", ".mat");
                return
            end
            obj = mlfourd.ImagingContext2(obj);
            fqfn = obj.fqfp + ".mat";
        end
        
        function populate_rawdata(opts)
            %% Populates CCIR_01211/rawdata from NIfTI previously generated by Nick Metcalf.
            %  N.B. alternative dates of FDG scanning.
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

        function time_align(fqfns, opts)
            %% aligns delay > 0 to delay == 0
            %  fqfns ~
            %  ["sub-108309_ses-20231204103722_trc-fdg_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz", ...
            %   "sub-108309_ses-20231204103722_trc-fdg_proc-delay300-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz"]
            %
            %  N.B.: for opts.noclobber, to avoid data corruption, aligning will not proceed if unaligned files
            %  exist, indicating aligning has already been done.

            arguments
                fqfns {mustBeText}
                opts.noclobber logical = true
                opts.use_static logical = true  % efficient if delay == 0 static is reliable
            end
            fqfns = convertCharsToStrings(fqfns);
            assert(all(arrayfun(@isfile, fqfns)))
            assert(all(arrayfun(@(x) contains(x, "-delay"), fqfns)))

            % no clobber
            unaligned = mglob(fullfile(unique(fileparts(fqfns)), "*-unaligned*"));
            if opts.noclobber && ~isempty(unaligned)
                return
            end

            % separate ref & its avgt
            idx_delay0 = find(contains(fqfns, "-delay0"));
            ref_ic = mlfourd.ImagingContext2(fqfns(idx_delay0));
            if opts.use_static
                ref_avgt_ic = mlfourd.ImagingContext2(strrep(ref_ic.fqfp, "createNiftiMovingAvgFrames", "createNiftiStatic"));
            else
                if ~isfile(ref_ic.fqfp + "_avgt.nii.gz")
                    ref_avgt_ic = ref_ic.timeAveraged();
                    ref_avgt_ic.save();
                else
                    ref_avgt_ic = mlfourd.ImagingContext2(ref_ic.fqfp + "_avgt.nii.gz");
                end
            end

            % rename non-ref files "-delay*" to "-unaligned*"
            fqfns1 = fqfns;
            fqfns1(idx_delay0) = [];
            unaligned_fqfns1 = strrep(fqfns1, "-delay", "-unaligned");
            json_fqfns1 = strrep(fqfns1, ".nii.gz", ".json");
            unaligned_json_fqfns1 = strrep(json_fqfns1, "-delay", "-unaligned");
            arrayfun(@(x, y) movefile(x, y), fqfns1, unaligned_fqfns1);
            arrayfun(@(x, y) movefile(x, y), json_fqfns1, unaligned_json_fqfns1);

            % construct non-ref avgt
            unaligned_avgt_fqfn = string();
            for idx = 1:length(unaligned_fqfns1)
                unaligned_ic__ = mlfourd.ImagingContext2(unaligned_fqfns1(idx));
                unaligned_avgt_ic__ = unaligned_ic__.timeAveraged();
                unaligned_avgt_ic__.save();
                unaligned_avgt_fqfn(idx) = unaligned_avgt_ic__.fqfn;
            end

            % flirt all non-ref avgt to ref avgt; apply transformations to time-series
            outs = strrep(unaligned_avgt_fqfn, ".nii.gz", "_on_ref.nii.gz");
            omats = strrep(unaligned_avgt_fqfn, ".nii.gz", "_on_ref.mat");
            for idx = 1:length(unaligned_avgt_fqfn)
                flirt = mlfsl.Flirt( ...
                    'in', unaligned_avgt_fqfn(idx), ...
                    'ref', ref_avgt_ic, ...
                    'out', outs(idx), ...
                    'omat', omats(idx), ...
                    'bins', 256, ...
                    'cost', 'mutualinfo', ...
                    'dof', 6, ...
                    'interp', 'trilinear', ...
                    'noclobber', false);
                if ~opts.noclobber || ~isfile(outs(idx))
                    % do expensive coreg.
                    flirt.flirt();
                    assert(isfile(outs(idx)))

                    % apply to time-series
                    flirt.in = unaligned_fqfns1(idx);
                    flirt.out = fqfns1(idx);
                    flirt.ref = ref_avgt_ic;
                    flirt.interp = 'trilinear';
                    flirt.applyXfm();
                    assert(isfile(fqfns1(idx)))
                end
            end

            % also align the static delay > 0
            mlvg.Lee2025.time_align_static(fqfns)
        end

        function time_align_static(fqfns)
            %% having completed time_align, also align static delay > 0 to static delay == 0
            %  fqfns ~
            %  ["sub-108309_ses-20231204103722_trc-fdg_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz", ...
            %   "sub-108309_ses-20231204103722_trc-fdg_proc-delay300-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz"]
            %
            %  N.B.: to avoid data corruption, aligning will not proceed if delayed_static & unaligned_delayed_static
            %  both exist, indicating aligning has already been done.

            arguments
                fqfns {mustBeText}
            end
            fqfns = convertCharsToStrings(fqfns);
            assert(all(arrayfun(@isfile, fqfns)))
            assert(all(arrayfun(@(x) contains(x, "-delay"), fqfns)))

            % find delay > 0
            all_matches = regexp(fqfns, "-delay(\d+)", "tokens", "once");  % returns cell array
            valid_matches = all_matches(~arrayfun(@isemptytext, all_matches));
            delay_values = cellfun(@(x) str2double(x{1}), valid_matches);
            delay_values = delay_values(delay_values ~= 0);
            delay_sec = delay_values(1);  % e.g., 30, 300

            % proceed only if there exists transformation:  delayed -> delay0
            fqfn_delay0 = fqfns(contains(fqfns, "-delay0"));
            fqfn_delayed = fqfns(contains(fqfns, "-delay" + delay_sec));
            prefix__ = strrep(fqfn_delayed, "-delay", "-unaligned");
            transform = myfileprefix(prefix__) + "_avgt_on_ref.mat";
            if ~isfile(transform)
                % e.g.: sub-108309_ses-20231204103722_trc-fdg_proc-unaligned300-BrainMoCo2-createNiftiMovingAvgFrames_avgt_on_ref.mat
                return
            end
            
            % apply transform (xfm) to static delay > 0
            delay0_static = strrep(fqfn_delay0, "createNiftiMovingAvgFrames", "createNiftiStatic");  % just needs to have correct image shape
            delayed_static = strrep(fqfn_delayed, "createNiftiMovingAvgFrames", "createNiftiStatic");
            unaligned_delayed_static = strrep(delayed_static, "-delay", "-unaligned");
            if isfile(delayed_static) && isfile(unaligned_delayed_static)
                % alignment was already done; don't clobber and confuse
                return
            end
            movefile(delayed_static, unaligned_delayed_static);
            delayed_static_j = strrep(delayed_static, ".nii.gz", ".json");
            unaligned_delayed_static_j = strrep(unaligned_delayed_static, ".nii.gz", ".json");
            copyfile(delayed_static_j, unaligned_delayed_static_j);
            flirt = mlfsl.Flirt( ...
                'in', unaligned_delayed_static, ...
                'ref', delay0_static, ...
                'out', delayed_static, ...
                'omat', transform, ...
                'bins', 256, ...
                'cost', 'mutualinfo', ...
                'dof', 6, ...
                'interp', 'trilinear', ...
                'noclobber', true);
            flirt.applyXfm();

            deleteExisting("*.log");
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
