classdef Lee2025 < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 04-Jul-2025 01:08:28 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.
   

    properties (Constant)
        % PARC_SCHAEF_TAG = "-ParcSchaeffer-reshape-to-schaeffer-schaeffer"
        % PARC_SCHAEF_TAG = "-ParcSchaeffer-invariant-schaeffer-schaeffer"
        PARC_SCHAEF_TAG = "-ParcSchaeffer-highsnr-schaeffer-schaeffer"
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
                opts.reference_tracer {mustBeTextScalar} = "fdg"
                opts.steps {mustBeNumericOrLogical} = 1
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
        
        function build_3dresample(fqfn, opts)
            %% builds 3dresample only on T1w

            arguments
                fqfn {mustBeFile}
                opts.noclobber logical = true
                opts.do_find_t1w logical = true
            end

            % fqfn is a file, but must be absolutely fully-qualified
            if ~contains(fqfn, pwd)
                fqfn = fullfile(pwd, fqfn);
            end

            % find T1w
            if opts.do_find_t1w && ~contains(mybasename(fqfn), "T1w")
                sub_path = extractBefore(fqfn, "/ses-");
                fqfns = mglob(fullfile(sub_path, "ses-*", "anat", "sub-*_ses-*_T1w_MPR_vNav_4e_RMS*.nii.gz"));
                if isempty(fqfns)
                    return
                end
                fqfn = fqfns(end);
            end

            % N.B. noclobber
            if opts.noclobber && contains(mybasename(fqfn), "orient-std")
                return
            end

            % directly build 3dresample without overhead of mlpipeline.ImagingMediator
            ic = mlfourd.ImagingContext2(fqfn);
            src_pth = fileparts(fqfn);
            derivs_pth = strrep(src_pth, "sourcedata", "derivatives");
            ensuredir(derivs_pth);

            pwd0 = pushd(src_pth);
            ic.afni_3dresample(orient_std=true);
            if ~strcmp(derivs_pth, ic.filepath)
                mysystem(sprintf("cp -f %s.* %s", ic.fqfileprefix, derivs_pth));
                mysystem(sprintf("rm -f %s.*", ic.fqfileprefix));
            end
            popd(pwd0);
        end

        function build_all(fqfn, opts)
            arguments
                fqfn {mustBeFile}
                opts.noclobber logical = true
            end

            % import mlsiemens.BrainMoCoBuilder
            % import mlvg.Lee2025.build_3dresample
            % import mlvg.Lee2025.repair_orient_std
            % import mlvg.Lee2025.repair_orient_std2
            % import mlvg.Lee2025.qc_orient_std
            % import mlvg.Lee2025.qc_schaeffer_parc
            % import mlvg.Lee2025Par.cluster_construct_pet_mipt
            % import mlvg.Lee2025Par.cluster_construct_pet_avgt
            % import mlvg.Lee2025Par.cluster_reflirt_t1w
            % import mlvg.Lee2025Par.cluster_time_align

            cd(fullfile(getenv("SINGULARITY_HOME"), "CCIR_01211"));
            build_workspace();  % in CCIR_01211/
            
            % N.B. mlsiemens.BrainMoCoBuilder methods for building derivatives/, sourcedata/
            build_3dresample(srcdata_t1w, noclobber=opts.noclobber);
            % ic = mlfourd.ImagingContext2();
            % ic = ic.blurred(3.5);
            % ic.save();

            % adjust objects in Schaefer2018_200Parcels_7Networks_order
            repair_orient_std(derivs_schaefer2018_reg);
            repair_orient_std2(derivs_schaefer2018_reg);
            qc_orient_std(derivs_t1w_os(4:end))
            T_qc_sch = qc_schaeffer_parc(derivs_schaefer2018_reg, do_repair=true);

            % intermediates for spatial normalization
            cluster_construct_pet_mipt(srcdata_fdg_delay0);
            cluster_construct_pet_avgt(srcdata_ho);
            cluster_reflirt_t1w("srcdata_ho.mat", globbing_var="srcdata_ho");
            cluster_reflirt_t1w("srcdata_fdg_delay0.mat", globbing_var="srcdata_fdg_delay0");
            cluster_reflirt_t1w("srcdata_oo_delay0.mat", globbing_var="srcdata_oo_delay0");
            cluster_reflirt_t1w("srcdata_co.mat", globbing_var="srcdata_co");

            % align delay0 with delay30|delay300
            cluster_time_align("srcdata_fdg_delay0.mat", globbing_var="srcdata_fdg_delay0");
            cluster_time_align("srcdata_oo_delay0.mat", globbing_var="srcdata_oo_delay0");

            % IDIFs
            build_mip_idif_finite(srcdata_oo_delay0)
            build_mip_idif_finite(srcdata_ho)
            build_mip_idif_finite(srcdata_co)
            build_mip_idif_finite(srcdata_fdg_delay0)

            % Schaefer and other parcels
            cluster_build_schaeffer_parc()
            build_schaeffer_delays()
            build_schaeffer_finite()

            % assemble Martin v1
            build_all_martin_v1_idif()

            % python models idif2025
            
            % update workspace for auditing
            build_workspace();  % in CCIR_01211/
        end

        function build_all_martin_v1_idif(mat_file)
            %% bulid v1 on Schaefer parcels using existing IDIFs for CO

            arguments
                mat_file {mustBeFile} = fullfile( ...
                    getenv("SINGULARITY_HOME"), "CCIR_01211", "srcdata_co.mat")
            end

            import mlkinetics.*

            ld = load(mat_file);
            globbed_co = ld.srcdata_co;

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

        function ic = build_mip_idif_finite(nmaf, opts)
            %% preferred method for building IDIF, superceding call_ifk

            arguments
                nmaf {mustBeText}
                opts.noclobber logical = false
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

        function ic = build_schaeffer_delays(nmaf, opts)
            %% Use with build_schaeffer_parc.
            %  Builds delayed NIfTI separately as needed.

            arguments
                nmaf {mustBeText}
                opts.noclobber logical = false
                opts.do_save logical = true
            end            

            ic = [];

            try
                % mark tracers with multiple delay sets
                if contains(nmaf, "trc-oo") || contains(nmaf, "trc-fdg")
                    nmaf = regexprep(nmaf, "-delay\d+", "-delay*");
                end

                % call base_case; save
                if isfile(nmaf)
                    ic = base_case(nmaf);
                    if opts.do_save
                        ic.save();
                    end
                    return
                end
                if contains(nmaf, "*")
                    nmafs = mglob(nmaf);
                    nmafs = nmafs(~contains(nmafs, "_avgt") & ~contains(nmafs, "_mipt"));
                    nmafs = natsort(nmafs);
                    ics = [];
                    for n = nmafs
                        ic = base_case(n);
                        if opts.do_save && ~isempty(ic)
                            ic.save();
                        end
                    end
                    return
                end
            catch ME
                handwarning(ME)
            end

            function ic = base_case(nmaf)
                ic = [];
                try

                    % check noclobber
                    [pth,fp] = myfileparts(nmaf);
                    pth_derivs = strrep(pth, "sourcedata", "derivatives");
                    fp_final = fp + mlvg.Lee2025.PARC_SCHAEF_TAG;
                    fqfn_final = fullfile(pth_derivs, fp_final + ".nii.gz");
                    if opts.noclobber && isfile(fqfn_final)
                        ic = mlfourd.ImagingContext2(fqfn_final);
                        return
                    end

                    % assign foundT1w
                    petMed = mlvg.Ccir1211Mediator.create(nmaf);
                    fp = lower(petMed.fileprefix);
                    foundT1w = mglob(fullfile(petMed.derivPetPath, sprintf("T1w_on_%s.nii.gz", fp)));
                    foundT1w = foundT1w(contains(foundT1w, extractBefore(mybasename(nmaf), "-BrainMoCo2")));
                    if isempty(foundT1w)
                        return
                    end

                    % assign foundRef
                    foundRef = mglob(fullfile(petMed.sourcePetPath, extractAfter(mybasename(foundT1w, withext=true), "T1w_on_")));
                    foundRef = foundRef(contains(foundRef, extractBefore(mybasename(nmaf), "-BrainMoCo2")));
                    if isempty(foundRef)
                        return
                    end
                    imagingReference = mlfourd.ImagingContext2(foundRef);
                    schaef_flirted_fqfn = strcat(petMed.fqfp, "-schaeffer.nii.gz");
                    schaef_flirted_fqfn = strrep(schaef_flirted_fqfn, "sourcedata", "derivatives");
                    target_fqfn = strrep(schaef_flirted_fqfn, "-schaeffer", mlvg.Lee2025.PARC_SCHAEF_TAG);
                    if isfile(target_fqfn)  % from ic1 = p.reshape_to_parc_fast(fqfn);
                        return
                    end

                    % flirt apply transform
                    omat = fullfile(petMed.derivPetPath, mybasename(foundT1w) + ".mat");
                    if ~isfile(omat)
                        % use alternative mat, sometimes created manually
                        % omat = mglob(fullfile(petMed.derivSubPath, "ses-*", "pet", sprintf("T1w_on_%s.mat", trc)));
                        omat = mglob(fullfile(petMed.derivSubPath, "ses-*", "pet", mybasename(foundT1w) + ".mat"));
                        omat = omat(1);
                    end
                    assert(isfile(omat))
                    flirt = mlfsl.Flirt( ...
                        'in', petMed.schaeffer_ic, ...
                        'ref', imagingReference, ...
                        'out', schaef_flirted_fqfn, ...
                        'omat', omat, ...
                        'bins', 4096, ...
                        'cost', 'normmi', ...
                        'interp', 'nearestneighbour', ...
                        'noclobber', true);
                    flirt.applyXfm();

                    % make bids kit, parc kit, parc, reshape to parc
                    bk = mlkinetics.BidsKit.create(bids_tags="ccir1211", bids_fqfn=schaef_flirted_fqfn);
                    pk = mlkinetics.ParcKit.create(bids_kit=bk, parc_tags="schaeffer-schaeffer");
                    p = pk.make_parc();
                    ic = p.reshape_to_parc_fast(nmaf);  % petMed.imagingContext, using mlvg.Lee2025.PARC_SCHAEF_TAG
                    ic.fqfn = fqfn_final;
                catch ME
                    handwarning(ME)
                end
            end
        end

        function ic = build_schaeffer_finite(nmaf, opts)
            %% Use with build_schaeffer_parc.
            %  Combines delays.
            %  Removes empty dynamic frames.

            arguments
                nmaf {mustBeText}
                opts.noclobber logical = false
                opts.do_save logical = true
            end

            ic = [];

            try
                % mark tracers with multiple delay sets
                if contains(nmaf, "trc-oo") || contains(nmaf, "trc-fdg")
                    nmaf = regexprep(nmaf, "-delay\d+", "-delay*");
                end

                % call base_case; save
                if isfile(nmaf)
                    ic = base_case(nmaf);
                    if opts.do_save
                        ic.save();
                    end
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
                    if opts.do_save
                        ic.save();
                    end
                    return
                end
            catch ME
                handwarning(ME)
            end

            function ic = base_case(nmaf)
                ic = [];
                try

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

                    % ensure finite
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
                    ic.fqfn = fqfn_final;
                catch ME
                    handwarning(ME);
                end
            end
        end

        function duration = build_schaeffer_parc(fqfn, opts)
            %% 1st:  build all delay\d+ 
            %  2nd:  build finite combining delays;
            %  e.g.,
            %  sub-108306_ses-20230227113853_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz ->
            %  sub-108306_ses-20230227113853_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames-ParcSchaeffer-reshape-to-schaeffer-schaeffer.nii.gz ->
            %  sub-108007_ses-20230227113853_trc-ho_proc-ParcSchaeffer-reshape-to-schaeffer-schaeffer-finite.nii.gz,
            %  ending with concatenation of delays when multiple delays exist.

            arguments
                fqfn {mustBeTextScalar}
                opts.do_plot logical = false
                opts.do_make_delays logical = true
                opts.do_make_finite logical = true
            end

            import mlkinetics.*

            duration = nan;
            try
                tic

                % 1st build schaeffer 2D NIfTI
                if opts.do_make_delays
                    mlvg.Lee2025.build_schaeffer_delays(fqfn, do_save=true);
                end

                % 2nd filter schaeffer 2D NIfTI to be finite
                if opts.do_make_finite
                    ic = mlvg.Lee2025.build_schaeffer_finite(fqfn, do_save=true);
                end

                if opts.do_plot
                    disp(ic.fqfn)
                    imagesc(ic)
                end

                duration = toc;
            catch ME
                handwarning(ME)
            end
        end

        function fqfn1 = construct_pet_avgt(fqfn, opts)
            %% work-around for failures of this.call_ifk()

            arguments
                fqfn {mustBeFile}
                opts.noclobber logical = true
            end

            % init
            [filepath,fileprefix,ext] = myfileparts(fqfn);
            filepath1 = strrep(filepath, "sourcedata", "derivatives");
            filename1 = fileprefix + "_avgt" + ext;
            fqfn1 = fullfile(filepath1, filename1);
            
            max_len_avgt = 60;

            % return existing fqfn1
            if opts.noclobber && isfile(fqfn1)
                return
            end

            % construct mipt
            ifc = mlfourd.ImagingFormatContext2(fqfn);
            img = ifc.img;
            L = min(size(img, 4), max_len_avgt);
            img = img(:,:,:,1:L);

            ifc.img = mean(img, 4, "omitmissing");
            ifc.fqfn = fqfn1;
            ifc.save();
        end

        function fqfn1 = construct_pet_mipt(fqfn, opts)
            %% work-around for failures of this.call_ifk()

            arguments
                fqfn {mustBeFile}
                opts.noclobber logical = true
            end

            % init
            [filepath,fileprefix,ext] = myfileparts(fqfn);
            filepath1 = strrep(filepath, "sourcedata", "derivatives");
            filename1 = fileprefix + "_mipt" + ext;
            fqfn1 = fullfile(filepath1, filename1);
            
            max_len_mipt = 30;
            minz_for_mip = 5;

            % return existing fqfn1
            if opts.noclobber && isfile(fqfn1)
                return
            end

            % construct mipt
            ifc = mlfourd.ImagingFormatContext2(fqfn);
            img = ifc.img;
            L = min(size(img, 4), max_len_mipt);
            img = img(:,:,:,1:L);
            img(:,:,1:minz_for_mip-1,:) = 0;
            img(:,:,end-minz_for_mip+1:end,:) = 0;

            ifc.img = max(img, [], 4);
            ifc.fqfn = fqfn1;
            ifc.save();
        end

        function fqfn1_final = deepmrseg_apply(fqfn, opts)
            arguments
                fqfn {mustBeFile}
                opts.blur {mustBeNumeric} = 7
            end

            % check fqfn is T1, check for existing DLICV
            assert(contains(mybasename(fqfn), "T1", IgnoreCase=true))            
            fqfn1 = strrep(fqfn, ".nii.gz", sprintf("_%s.nii.gz", mlvg.Ccir1211Bids.DLICV_TAG));
            if ~isempty(opts.blur)
                fqfn1_final = strrep(fqfn1, ".nii.gz", sprintf("_b%g.nii.gz", round(10*max(opts.blur))));
            else
                fqfn1_final = fqfn1;
            end
            if isfile(fqfn1_final)
                return
            end

            % call deepmrseg_apply
            ic1 = mlpipeline.Bids.deepmrseg_apply(fqfn, fqfn1);
            if ~isempty(opts.blur)
                ic1.save();
                ic1 = ic1.blurred(opts.blur);
            end
            ic1.fqfn = fqfn1_final;
            ic1.save();
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
                opts.use_avgt logical = true
            end

            subpth = extractBefore(fileparts(other_nii), filesep + "ses-");

            if opts.use_avgt
                subpth = strrep(subpth, "sourcedata", "derivatives");
                globbed = mglob(fullfile(subpth, "ses-*", "pet", ...
                    "sub-*_ses-*_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_avgt.nii.gz"));
                assert(~isempty(globbed))
                assert(isscalar(globbed(opts.selection)))
                fqfn = globbed(opts.selection);
                return
            end

            globbed = mglob(fullfile(subpth, "ses-*", "pet", ...
                "sub-*_ses-*_trc-ho_proc-delay0-BrainMoCo2-createNiftiStatic.nii.gz"));
            assert(~isempty(globbed))
            assert(isscalar(globbed(opts.selection)))
            fqfn = globbed(opts.selection);
        end

        function fqfn = find_t1w(other_nii, opts)
            %% e.g., sub-108333_ses-20241122094951_T1w_MPR_vNav_4e_RMS_orient-std.nii.gz

            arguments
                other_nii {mustBeFile}
                opts.tag {mustBeTextScalar} = "_b35"
            end

            subpth = extractBefore(fileparts(other_nii), filesep + "ses-");
            subpth = strrep(subpth, "sourcedata", "derivatives");
            globbed = mglob(fullfile(subpth, "ses-*", "anat", ...
                sprintf("sub-*_ses-*_T1w_MPR_vNav_4e_RMS*_orient-std%s.nii.gz", opts.tag)));
            assert(~isempty(globbed))
            globbed = globbed(1);
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
                opts.use_avgt logical = true
            end

            subpth = extractBefore(fileparts(other_nii), filesep + "ses-");
            subpth = strrep(subpth, "sourcedata", "derivatives");
            if opts.use_avgt
                globbed = mglob(fullfile(subpth, "ses-*", "pet", ...
                    "T1w_on_*_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames_avgt.nii.gz"));
            else
                globbed = mglob(fullfile(subpth, "ses-*", "pet", ...
                    "T1w_on_*_trc-ho_proc-delay0-BrainMoCo2-createNiftiStatic.nii.gz"));            
            end
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
            arguments
                fqfn {mustBeFile}
                opts.specialize_for_tracer logical = true
                opts.noclobber logical = true
                opts.cost {mustBeTextScalar} = "normmi"
                opts.use_dlicv logical = false
            end

            if contains(fqfn, "sub-108259_ses-20230731133714")
                % OO administered but scanning performed with CO protocol
                opts.specialize_for_tracer = false;
            end
            if opts.specialize_for_tracer
                if contains(fqfn, "trc-ho")
                    mlvg.Lee2025.flirt_t1w_on_ho(fqfn, noclobber=opts.noclobber, cost=opts.cost);
                    return
                end
                if contains(fqfn, "trc-co")
                    mlvg.Lee2025.flirt_t1w_on_co(fqfn, noclobber=opts.noclobber, cost=opts.cost);
                    return
                end
                if contains(fqfn, "trc-oo")
                    mlvg.Lee2025.flirt_t1w_on_oo(fqfn, noclobber=opts.noclobber, cost=opts.cost);
                    return
                end
            end

            import mlvg.Lee2025.mat

            fqfn = strrep(fqfn, "createNiftiMovingAvgFrames", "createNiftiStatic");

            % flirt fdg_mipt -> co_mipt
            t1w_fqfn = mlvg.Lee2025.find_t1w(fqfn);
            [pth,fp] = myfileparts(fqfn);
            dlicv_fqfn = [];
            if opts.use_dlicv
                dlicv_glob = mglob(fullfile(fileparts(t1w_fqfn), "*_DLICV_b70.nii.gz"));
                if ~isempty(dlicv_glob)
                    dlicv_fqfn = dlicv_glob(1);
                end
            end
            pth_derivs = strrep(pth, "sourcedata", "derivatives");
            t1w_on_tracer = fullfile(pth_derivs, "T1w_on_" + fp + ".nii.gz");
            flirt = mlfsl.Flirt( ...
                'in', t1w_fqfn, ...
                'inweight', dlicv_fqfn, ...
                'ref', fqfn, ...
                'out', t1w_on_tracer, ...
                'omat', mat(t1w_on_tracer), ...
                'bins', 4096, ...
                'cost', opts.cost, ...
                'searchrx', 20, ...
                'searchry', 20, ...
                'searchrz', 20, ...
                'dof', 6, ...
                'interp', 'spline', ...
                'noclobber', false);
            ensuredir(fileparts(t1w_on_tracer))
            if ~opts.noclobber || ~isfile(t1w_on_tracer)
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
                opts.cost {mustBeTextScalar} = "normmi"
            end

            co_fqfn = strrep(co_fqfn, "createNiftiMovingAvgFrames", "createNiftiStatic");
            co_fqfn = strrep(co_fqfn, "consoleDynamic", "consoleStatic");

            import mlvg.Lee2025.find_fdg_mipt
            import mlvg.Lee2025.find_t1w
            import mlvg.Lee2025.find_t1w_on_fdg
            import mlvg.Lee2025.mat

            fdg_mipt = find_fdg_mipt(co_fqfn);
            fdg_mipt_on_co = myfileprefix(fdg_mipt) + "_on_" + mybasename(co_fqfn) + ".nii.gz";

            % flirt fdg_mipt -> co_mipt
            flirt = mlfsl.Flirt( ...
                'in', fdg_mipt, ...
                'ref', co_fqfn, ...
                'out', fdg_mipt_on_co, ...
                'omat', mat(fdg_mipt_on_co), ...
                'bins', 4096, ...
                'cost', opts.cost, ...
                'searchrx', 20, ...
                'searchry', 20, ...
                'searchrz', 20, ...
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

        function flirt_t1w_on_ho(ho_fqfn, opts)
            %% intended for failures of MipIdif;
            %  try again using ho_static -> oo_delay30_avgt, then applyXfm(T1w);
            %  e.g., sub-108034_ses-20230717133405_trc-oo_proc-delay30-BrainMoCo2-createNiftiStatic.nii.gz

            arguments
                ho_fqfn {mustBeFile}
                opts.noclobber logical = true
                opts.cost {mustBeTextScalar} = "normmi"
                opts.use_dlicv logical = false
            end

            import mlvg.Lee2025.mat

            ho_fqfn = strrep(ho_fqfn, "sourcedata", "derivatives");
            ho_fqfn = strrep(ho_fqfn, "createNiftiMovingAvgFrames", "createNiftiMovingAvgFrames_avgt");

            % flirt fdg_mipt -> co_mipt
            t1w_fqfn = mlvg.Lee2025.find_t1w(ho_fqfn);
            [pth,fp] = myfileparts(ho_fqfn);
            dlicv_fqfn = [];
            if opts.use_dlicv
                dlicv_glob = mglob(fullfile(fileparts(t1w_fqfn), "*_DLICV_b70.nii.gz"));
                if ~isempty(dlicv_glob)
                    dlicv_fqfn = dlicv_glob(1);
                end
            end
            pth_derivs = strrep(pth, "sourcedata", "derivatives");
            t1w_on_tracer = fullfile(pth_derivs, "T1w_on_" + fp + ".nii.gz");
            flirt = mlfsl.Flirt( ...
                'in', t1w_fqfn, ...
                'inweight', dlicv_fqfn, ...
                'ref', ho_fqfn, ...
                'out', t1w_on_tracer, ...
                'omat', mat(t1w_on_tracer), ...
                'bins', 4096, ...
                'cost', opts.cost, ...
                'searchrx', 20, ...
                'searchry', 20, ...
                'searchrz', 20, ...
                'dof', 6, ...
                'interp', 'spline', ...
                'noclobber', false);
            ensuredir(fileparts(t1w_on_tracer))
            if ~opts.noclobber || ~isfile(t1w_on_tracer)
                flirt.flirt();
                assert(isfile(t1w_on_tracer))
            end
        end

        function flirt_t1w_on_oo(oo_fqfn, opts)
            %% intended for failures of MipIdif;
            %  try again using ho_static -> oo_delay30_avgt, then applyXfm(T1w);
            %  e.g., sub-108034_ses-20230717133405_trc-oo_proc-delay30-BrainMoCo2-createNiftiStatic.nii.gz

            arguments
                oo_fqfn {mustBeFile}
                opts.noclobber logical = true
                opts.cost {mustBeTextScalar} = "normmi"
            end

            import mlvg.Lee2025.find_ho_avgt
            import mlvg.Lee2025.find_t1w
            import mlvg.Lee2025.find_t1w_on_ho
            import mlvg.Lee2025.mat

            oo_fqfn = strrep(oo_fqfn, "createNiftiMovingAvgFrames", "createNiftiStatic");

            try
                ho_avgt = mlvg.Lee2025.find_ho_avgt(oo_fqfn);
                ho_avgt_on_oo = myfileprefix(ho_avgt) + "_on_" + mybasename(oo_fqfn) + ".nii.gz";
                ho_avgt_on_oo = strrep(ho_avgt_on_oo, "sourcedata", "derivatives");
            catch ME
                fprintf("mlvg:RuntimeError: %s %s", stackstr(), ME.message);
                mlvg.Lee2025.flirt_t1w(oo_fqfn, noclobber=opts.noclobber, cost=opts.cost, ...
                    specialize_for_tracer=false);
                return
            end

            % flirt ho_avgt -> oo
            flirt = mlfsl.Flirt( ...
                'in', ho_avgt, ...
                'ref', oo_fqfn, ...
                'out', ho_avgt_on_oo, ...
                'omat', mat(ho_avgt_on_oo), ...
                'bins', 4096, ...
                'cost', opts.cost, ...
                'searchrx', 20, ...
                'searchry', 20, ...
                'searchrz', 20, ...
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

        function resultTable = parse_flirt_results(jsonFilenames)
            % PARSE_FLIRT_RESULTS Parse FSL flirt co-registration JSON files into a table
            %
            % Syntax:
            %   resultTable = parse_flirt_results(jsonFilenames)
            %
            % Input:
            %   jsonFilenames - String array of JSON filenames
            %                   Can be a single string, char array, or string array
            %                   Can be relative or absolute paths
            %
            % Output:
            %   resultTable - MATLAB table with columns:
            %                 - json_file: basename of the JSON file
            %                 - transformation_file: basename of the transformation matrix file
            %                 - final_cost: final cost value from co-registration
            %
            % Example:
            %   % Single file
            %   T = parse_flirt_results("result.json");
            %
            %   % Multiple files using string array
            %   files = ["result1.json", "/path/to/result2.json", "data/result3.json"];
            %   T = parse_flirt_results(files);
            %
            %   % Using string constructor
            %   files = string({'result1.json', 'result2.json', 'result3.json'});
            %   T = parse_flirt_results(files);
            %
            % Author: Claude Opus 4.1
            % Date: 2025/09/28 11:34

            % Convert input to string array if necessary
            if ischar(jsonFilenames)
                jsonFilenames = string(jsonFilenames);
            elseif iscell(jsonFilenames)
                jsonFilenames = string(jsonFilenames);
            elseif ~isstring(jsonFilenames)
                error('Input must be a string array, char array, or cell array of filenames');
            end

            % Ensure column vector for consistency
            jsonFilenames = jsonFilenames(:);

            % Initialize string arrays and numeric array to store results
            nFiles = length(jsonFilenames);
            json_files = strings(nFiles, 1);
            transformation_files = strings(nFiles, 1);
            final_costs = nan(nFiles, 1);
            final_sum_abs = nan(nFiles, 1);

            % Process each JSON file
            for i = 1:nFiles
                filename = jsonFilenames(i);

                % Check if file exists
                if ~isfile(filename)
                    warning('File not found: %s. Skipping...', filename);
                    json_files(i) = mybasename(filename);
                    transformation_files(i) = "N/A";
                    final_costs(i) = NaN;
                    final_sum_abs(i) = NaN;
                    continue;
                end

                try
                    % Read and parse JSON file
                    jsonText = fileread(filename);
                    jsonData = jsondecode(jsonText);

                    % Store the basename of the JSON file
                    json_files(i) = mybasename(filename);

                    % Extract transformation file and final cost
                    % Handle potential variations in JSON structure
                    if isfield(jsonData, 'mlfsl_Flirt')
                        flirtData = jsonData.mlfsl_Flirt;

                        % Extract transformation file
                        if isfield(flirtData, 'cost_final') && ...
                                isfield(flirtData.cost_final, 'init')
                            transformation_files(i) = mybasename(string(flirtData.cost_final.init));
                        else
                            warning('Transformation file (init) not found in %s', filename);
                            transformation_files(i) = "N/A";
                        end

                        % Extract final cost
                        if isfield(flirtData, 'cost_final') && ...
                                isfield(flirtData.cost_final, 'cost')
                            final_costs(i) = flirtData.cost_final.cost;
                        else
                            warning('Final cost value not found in %s', filename);
                            final_costs(i) = NaN;
                        end
                    else
                        warning('mlfsl_Flirt field not found in %s', filename);
                        transformation_files(i) = "N/A";
                        final_costs(i) = NaN;
                    end

                    % dipsum(abs(ImagingContext2)) identifies spatial normalizations 
                    % placing brain outside the FOV
                    filename_niigz = strrep(filename, ".json", ".nii.gz");
                    ic = mlfourd.ImagingContext2(filename_niigz);
                    final_sum_abs(i) = dipsum(abs(ic));

                catch ME
                    warning('Error processing file %s: %s', filename, ME.message);
                    json_files(i) = mybasename(filename);
                    transformation_files(i) = "N/A";
                    final_costs(i) = NaN;
                    final_sum_abs(i) = NaN;
                end
            end

            % Create the output table
            resultTable = table(json_files, transformation_files, final_costs, final_sum_abs, ...
                'VariableNames', {'json_file', 'transformation_file', 'final_cost', 'final_sum_voxel'});

            % Display summary
            fprintf('\n=== Processing Summary ===\n');
            fprintf('Total files processed: %d\n', nFiles);
            fprintf('Successfully parsed: %d\n', sum(~isnan(final_costs)));
            fprintf('Failed/Missing: %d\n', sum(isnan(final_costs)));

            % Display warning about NaN values if any
            if any(isnan(final_costs))
                fprintf('\nNote: NaN values in final_cost indicate parsing errors or missing data.\n');
                fprintf('Check warnings above for details.\n');
            end
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

        function qc_orient_std(fqfns)
            arguments
                fqfns {mustBeText}
            end

            for f = fqfns
                ic = mlfourd.ImagingContext2(f);
                f1 = mglob(fullfile( ...
                    extractBefore(f, "/ses-"), "ses-*", "**", ...
                    "Schaefer2018_200Parcels_7Networks_order_T1_registered.nii.gz"));
                for f2 = f1
                    fprintf("%s: %s, %s\n", stackstr(), f, f2);
                    ic.view_qc(f2);
                end
            end
        end

        function T = qc_schaeffer_parc(fqfns, opts)
            arguments
                fqfns {mustBeText}
                opts.fqfn_ref {mustBeTextScalar} = ""
                opts.do_repair logical = false  % sets 80 -> 0; uses proxies for missing 30, 62
            end
            if isemptytext(opts.fqfn_ref)
                opts.fqfn_ref = fqfns(1);
            end

            fqfns = ascol(fqfns);
            num_parcs = nan(size(fqfns));
            xor = cell(size(fqfns));
            ifc_ref = mlfourd.ImagingFormatContext2(opts.fqfn_ref);
            parcs_ref = unique(ifc_ref.img);
            for idx = 1:length(fqfns)
                ifc = mlfourd.ImagingFormatContext2(fqfns(idx));
                if opts.do_repair
                    ifc = repair_parcs(ifc);
                end
                parcs = unique(ifc.img);
                num_parcs(idx) = numel(setxor(parcs, 0));
                xor{idx} = setxor(parcs, parcs_ref);
            end

            T = table(fqfns, num_parcs, xor);

            function ifc = repair_parcs(ifc)
                mandatory = [30, 62];  % Left-vessel, Right-vessel
                proxy = [31, 63];  % L/R-choroid-plexus 
                expendable = 80;  % non-WM-hypointensities
                do_save = false;
                ifc_bak = copy(ifc);

                % recover mandatory from sampling of proxy
                for midx = 1:length(mandatory)
                    m = mandatory(midx);
                    if sum(ifc.img == m, "all") == 0
                        select = ifc.img == proxy(midx);

                        % Find true indices, randomly select half, keep only selected trues
                        trueIdx = find(select);
                        keepIdx = trueIdx(randperm(length(trueIdx), floor(length(trueIdx)/2)));
                        select_ = false(size(select));
                        select_(keepIdx) = true;

                        ifc.img(select_) = m;
                        do_save = true;
                    end
                end

                % annihilate expendable
                for e = expendable
                    if sum(ifc.img == e, "all") > 0
                        ifc.img(ifc.img == e) = 0;
                        do_save = true;
                    end
                end
                if do_save
                    fprintf("%s: repairing %s\n", stackstr(), ifc_bak.fqfn);
                    ifc_bak.selectImagingFormatTool();
                    ifc_bak.fileprefix = ifc_bak.fileprefix + "_backup";
                    ifc_bak.save();  % overwrites existing backup
                    ifc.save();  % overwrites
                end
            end
        end

        function repair_orient_std(t1w_os)
            arguments
                t1w_os {mustBeText}
            end

            for t = asrow(t1w_os)
                if ~endsWith(t, "orient-std.nii.gz")
                    continue
                end
                if ~contains(t, pwd)
                    t1w_fqfn = fullfile(pwd, t);
                else
                    t1w_fqfn = t;
                end
                t1w_ifc = mlfourd.ImagingFormatContext2(t1w_fqfn);
                registered = copy(t1w_ifc);

                sch_oss = mglob(fullfile( ...
                    extractBefore(t1w_fqfn, "/ses-"), "ses-*", "**", ...
                    "Schaefer2018_200Parcels_7Networks_order_T1_orient-std.nii.gz"));
                for s = sch_oss
                    sch_os = mlfourd.ImagingFormatContext2(s);
                    registered.img = sch_os.img;
                    registered.filepath = sch_os.filepath;
                    registered.fileprefix = strrep(sch_os.fileprefix, "orient-std", "registered");
                    registered.save();
                end
                
            end
        end

        function repair_orient_std2(sch_reg, opts)
            arguments
                sch_reg {mustBeFile}
                opts.noclobber logical = false
                opts.use_orig_complete logical = false
            end

            t1w_os = mlvg.Lee2025.find_t1w(sch_reg, tag="");
            assert(isfile(t1w_os))
            orig = fullfile(fileparts(sch_reg), "orig.nii.gz");
            assert(isfile(orig))
            if opts.use_orig_complete
                sch_orig = strrep(sch_reg, "T1_registered", "orig_complete");
            else
                sch_orig = strrep(sch_reg, "T1_registered", "orig");
            end
            assert(contains(sch_orig, "orig"))
            assert(isfile(sch_orig))
            orig_on_t1w = fullfile(fileparts(sch_reg), sprintf("orig_on_%s.nii.gz", mybasename(t1w_os)));
            orig_on_t1w_mat = mlvg.Lee2025.mat(orig_on_t1w);

            % flirt orig -> T1w
            flirt = mlfsl.Flirt( ...
                'in', orig, ...
                'ref', t1w_os, ...
                'out', orig_on_t1w, ...
                'omat', orig_on_t1w_mat, ...
                'bins', 2048, ...
                'cost', 'corratio', ...
                'dof', 6, ...
                'interp', 'spline', ...
                'noclobber', false);
            if ~opts.noclobber || ~isfile(orig_on_t1w)
                flirt.flirt();
            end
            assert(isfile(orig_on_t1w))

            % apply xfm:  sch_reg -> sch_reg
            sch_reg_bak = strrep(sch_reg, ".nii", "_bak.nii");
            copyfile_niigz(sch_reg, sch_reg_bak);
            assert(isfile(sch_reg_bak))
            flirt.in = sch_orig;
            flirt.out = sch_reg;
            flirt.interp = "nearestneighbour";
            flirt.applyXfm();
            assert(isfile(sch_reg))
        end

        function fqfns1 = rename_delay(fqfns, opts)
            %% mv *delay*{nii.gz,json} to *unaligned*{nii.gz,json}

            arguments
                fqfns {mustBeText}
                opts.noclobber logical = true
                opts.do_revert logical = false  % mv *unaligned* to *delay*
            end

            % safeguard NIfTIs with delay0
            fqfns = fqfns(~contains(fqfns, "delay0"));

            if ~opts.do_revert
                assert(all(contains(fqfns, "-delay")))
                fqfns1 = string();
                for fidx = 1:length(fqfns)
                    fqfn = fqfns(fidx);
                    if ~isfile(fqfn)
                        continue
                    end
                    unaligned_fqfn = strrep(fqfn, "-delay", "-unaligned");
                    if opts.noclobber && isfile(unaligned_fqfn)
                        continue
                    end
                    arrayfun(@(x, y) movefile_niigz(x, y), fqfn, unaligned_fqfn, UniformOutput=false);
                    fqfns1(fidx) = unaligned_fqfn;
                end
                return
            end
            
            % reversion is more complex
            % -- fqfns must be all unaligned
            % -- each of fqfns must exist as a file
            % -- delete unaligned tagged "_avgt", including ".mat"
            % -- delete delayed that have corresponding unaligned
            % -- mv unaligned to delayed, to conclude reversion

            assert(all(contains(fqfns, "-unaligned")))
            fqfns1 = string();
            for fidx = 1:length(fqfns)
                fqfn = fqfns(fidx);
                if ~isfile(fqfn)
                    continue
                end
                if contains(fqfn, "_avgt")
                    delete_niigz(fqfn, exts=[".nii.gz", ".json", ".log", ".mat"]);
                    continue
                end
                delayed_fqfn = strrep(fqfn, "-unaligned", "-delay");
                if isfile(delayed_fqfn)
                    delete_niigz(delayed_fqfn);
                end
                arrayfun(@(x, y) movefile_niigz(x, y), fqfn, delayed_fqfn, UniformOutput=false);
                fqfns1(fidx) = delayed_fqfn;
            end
        end

        function [xor,re,re2] = setxor_sub_ses(filearr, filearr2)
            %% set xor matching only sub & ses

            arguments
                filearr {mustBeText}
                filearr2 {mustBeText}
            end

            re = regexp(filearr, "\S+/(?<subses>sub-\d{6}/ses-(\d{8}|\d{14}))/\S+", "names");
            re2 = regexp(filearr2, "\S+/(?<subses>sub-\d{6}/ses-(\d{8}|\d{14}))/\S+", "names");
            ssarr = string(cellfun(@(x) x.subses, re(~isempty(re)), UniformOutput=false));
            ssarr2 = string(cellfun(@(x) x.subses, re2(~isempty(re2)), UniformOutput=false));
            xor = setxor(ssarr, ssarr2);
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
                opts.use_static logical = true  % efficient if delay == 0 static is reliable as reference
                opts.use_refweight logical = false
            end
            fqfns = convertCharsToStrings(fqfns);
            assert(all(arrayfun(@isfile, fqfns)))
            assert(all(arrayfun(@(x) contains(x, "-delay"), fqfns)))

            import mlvg.Lee2025.rename_delay;

            % no clobber for existing files "unaligned"
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

            % rename non-ref files "-delay*" to "-unaligned*"; copy *.json
            fqfns1 = fqfns;
            fqfns1(idx_delay0) = [];
            unaligned_fqfns1 = rename_delay(fqfns1);
            
            % construct non-ref avgt
            unaligned_avgt_fqfn = string();
            for idx = 1:length(unaligned_fqfns1)
                % *delay*_avgt.nii.gz may exist
                avgt_fqfn = strrep(fqfns1(idx), ".nii.gz", "_avgt.nii.gz");
                if isfile(avgt_fqfn)
                    % avoid expensive construction of _avgt
                    unaligned_avgt_fqfn(idx) = rename_delay(avgt_fqfn);
                    continue
                end

                % create de novo *unaligned*_avgt.nii.gz
                unaligned_ic__ = mlfourd.ImagingContext2(unaligned_fqfns1(idx));
                unaligned_avgt_ic__ = unaligned_ic__.timeAveraged();
                unaligned_avgt_ic__.save();
                unaligned_avgt_fqfn(idx) = unaligned_avgt_ic__.fqfn;
            end

            % consider using refweight
            if opts.use_refweight
                try
                    % gather info
                    deriv_pet_path = strrep(ref_avgt_ic.filepath, "sourcedata", "derivatives");
                    deriv_sub_path = extractBefore(deriv_pet_path, "/ses-");
                    deriv_anat_path = mglob(fullfile(deriv_sub_path, "ses-*", "anat"));
                    deriv_anat_path = deriv_anat_path(1);
                    weight = mglob(fullfile(deriv_anat_path, "sub-*_ses-*_T1w*_DLICV_b*.nii.gz"));
                    weight = weight(end);
                    t1w_on_pet_fp = mglob(fullfile(deriv_pet_path, "T1w_on_sub-*_ses-*_trc*delay30*.nii.gz"));
                    t1w_on_pet_fp = mybasename(t1w_on_pet_fp(end));
                    re = regexp(mybasename(weight), "(?<dlicv>_DLICV_b\d+)$", "names");
                    refweight = fullfile(deriv_pet_path, t1w_on_pet_fp + re.dlicv + ".nii.gz");
                    mat = mglob(fullfile(deriv_sub_path, "ses-*", "pet", t1w_on_pet_fp + ".mat"));
                    assert(isfile(mat))

                    % applyXfm
                    flirt = mlfsl.Flirt( ...
                        'in', weight, ...
                        'ref', ref_avgt_ic, ...
                        'out', refweight, ...
                        'omat', mat, ...
                        'bins', 4096, ...
                        'cost', 'mutualinfo', ...
                        'searchrx', 10, ...
                        'searchry', 10, ...
                        'searchrz', 10, ...
                        'dof', 6, ...
                        'interp', 'trilinear', ...
                        'noclobber', false);
                    flirt.applyXfm();
                catch ME
                    handwarning(ME)
                    refweight = [];
                end
            else
                refweight = [];
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
                    'bins', 4096, ...
                    'cost', 'mutualinfo', ...
                    'searchrx', 10, ...
                    'searchry', 10, ...
                    'searchrz', 10, ...
                    'dof', 6, ...
                    'refweight', refweight, ...
                    'interp', 'spline', ...
                    'noclobber', true);
                if ~opts.noclobber || ~isfile(outs(idx))
                    % do expensive coreg.
                    flirt.flirt();
                    assert(isfile(outs(idx)))

                    % apply to time-series
                    flirt.in = unaligned_fqfns1(idx);
                    flirt.out = fqfns1(idx);
                    flirt.ref = ref_avgt_ic;
                    flirt.interp = 'spline';
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
            movefile_niigz(delayed_static, unaligned_delayed_static);
            %delayed_static_j = strrep(delayed_static, ".nii.gz", ".json");
            %unaligned_delayed_static_j = strrep(unaligned_delayed_static, ".nii.gz", ".json");
            %copyfile_niigz(delayed_static_j, unaligned_delayed_static_j);
            flirt = mlfsl.Flirt( ...
                'in', unaligned_delayed_static, ...
                'ref', delay0_static, ...
                'out', delayed_static, ...
                'omat', transform, ...
                'bins', 4096, ...
                'cost', 'mutualinfo', ...
                'searchrx', 10, ...
                'searchry', 10, ...
                'searchrz', 10, ...
                'dof', 6, ...
                'interp', 'spline', ...
                'noclobber', true);
            flirt.applyXfm();

            deleteExisting(fullfile(fileparts(delayed_static), "*.log"));
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
