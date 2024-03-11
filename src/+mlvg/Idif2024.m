classdef Idif2024 < handle & mlsystem.IHandle
    %% Provides simpler, high-level, re-usable application programming interfaces for software supporting
    %  Lee, et al., to be submitted to Physics in Biology & Medicine, 2024.
    %
    %  See also:  mlvg.Lee2024, PycharmProjects/dynesty/idif2024, Singularity/CCIR_01222.
    %  
    %  Created 06-Mar-2024 23:13:15 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 23.2.0.2515942 (R2023b) Update 7 for MACA64.  Copyright 2024 John J. Lee.
    

    properties (Constant)
        LC = 0.81  % lumped constant
        RECOVERY_COEFFICIENT = 1.8509  % determined from CO for Biograph Vision 600
        TIME_APPEND_SUFFIXES = ["_timeAppend-165", "_timeAppend-80", "_timeAppend-4", ""]
    end

    properties (Dependent)
        derivatives_path
        glc  % mg/dL
        o2_content  % containers.map:  sub-id =: o2 content ~ [0, 1] in mL/mL
        sourcedata_path
    end

    methods  %% GET
        function g = get.derivatives_path(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "CCIR_01211", "derivatives");
        end
        function g = get.glc(this)
            if ~isempty(this.glc_)
                g = this.glc_;
                return
            end

            this.glc_ = containers.Map();
            this.glc_("sub-108293") = 114;  % mg/dL
            this.glc_("sub-108237") = 105;
            this.glc_("sub-108254") = 99;
            this.glc_("sub-108250") = 89;
            this.glc_("sub-108284") = 109;
            this.glc_("sub-108306") = 107;
            g = this.glc_;
        end
        function g = get.o2_content(this)
            if ~isempty(this.o2_content_)
                g = this.o2_content_;
                return
            end

            this.o2_content_ = containers.Map();
            this.o2_content_("sub-108293") = 0.179;  % mL/mL
            this.o2_content_("sub-108237") = 0.184;
            this.o2_content_("sub-108254") = 0.159;
            this.o2_content_("sub-108250") = 0.189;
            this.o2_content_("sub-108284") = 0.184;
            this.o2_content_("sub-108306") = 0.179;
            g = this.o2_content_;
        end
        function g = get.sourcedata_path(this)
            g = strrep(this.derivatives_path, "derivatives", "sourcedata");
        end
    end

    methods

        %% physiological data objects

        function obj = cbv(this, opts)
            %% mL/mL median over all available data

            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            ic = this.martinv1(input_func=opts.input_func, typeclass="nifti", stats=opts.stats);
            ic.fileprefix = strrep(ic.fileprefix, "martinv1", "cbv");
            obj = this.as(ic, opts.typeclass);
        end

        function obj = cbf(this, opts)
            %% mL/mL/min median over all available data

            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            ic = this.raichleks(input_func=opts.input_func, typeclass="nifti", stats=opts.stats);
            ifc = ic.imagingFormat;
            ifc.img = ifc.img(:,1) .* 60;
            ifc.fileprefix = strrep(ifc.fileprefix, "raichleks", "cbf");
            ic = mlfourd.ImagingContext2(ifc);
            obj = this.as(ic, opts.typeclass);
        end

        function obj = lambda(this, opts)
            %% partition coefficient for water, mL/mL

            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            ic = this.raichleks(input_func=opts.input_func, typeclass="nifti", stats=opts.stats);
            ifc = ic.imagingFormat;
            ifc.img = ifc.img(:,2);
            ifc.fileprefix = strrep(ifc.fileprefix, "raichleks", "lambda");
            ic = mlfourd.ImagingContext2(ifc);
            obj = this.as(ic, opts.typeclass);
        end

        function obj = ps(this, opts)
            %% permeability-surface area product, mL/mL/min

            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            ic = this.raichleks(input_func=opts.input_func, typeclass="nifti", stats=opts.stats);
            ifc = ic.imagingFormat;
            ifc.img = ifc.img(:,3) .* 60;
            ifc.fileprefix = strrep(ifc.fileprefix, "raichleks", "ps");
            ic = mlfourd.ImagingContext2(ifc);
            obj = this.as(ic, opts.typeclass);
        end

        function obj = oef(this, opts)
            %% mL/mL

            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            ic = this.mintunks(input_func=opts.input_func, typeclass="nifti", stats=opts.stats);
            ifc = ic.imagingFormat;
            ifc.img = ifc.img(:,1);
            ifc.fileprefix = strrep(ifc.fileprefix, "mintunks", "oef");
            ic = mlfourd.ImagingContext2(ifc);
            obj = this.as(ic, opts.typeclass);
        end

        function obj = vcapillary(this, opts)
            %% r"$v_{post} + 0.5 v_{cap}$", mL/mL

            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            ic = this.mintunks(input_func=opts.input_func, typeclass="nifti", stats=opts.stats);
            ifc = ic.imagingFormat;
            ifc.img = ifc.img(:,2);
            ifc.fileprefix = strrep(ifc.fileprefix, "mintunks", "vcapillary");
            ic = mlfourd.ImagingContext2(ifc);
            obj = this.as(ic, opts.typeclass);
        end

        function obj = fwatermetab(this, opts)
            %% r"frac. water of metab.", mL/mL

            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            ic = this.mintunks(input_func=opts.input_func, typeclass="nifti", stats=opts.stats);
            ifc = ic.imagingFormat;
            ifc.img = ifc.img(:,3);
            ifc.fileprefix = strrep(ifc.fileprefix, "mintunks", "fwatermetab");
            ic = mlfourd.ImagingContext2(ifc);
            obj = this.as(ic, opts.typeclass);
        end

        function [obj,mg] = cmro2(this, opts)
            %% mmol / L / min
            %  N.B.:  cmro2 := cmro2 * 44.64; mmol O_2 / min / L := mL O_2 / min / hg

            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            % for each available ic_mintunks:
            %     aufbau cmro2 from this.o2_content, ifc_mintunks, ic_raichleks
            [ic_mintunks,mg] = this.mintunks(input_func=opts.input_func, typeclass="nifti", stats="");
            ifc_mintunks = ic_mintunks.imagingFormat;  % img ~ Nvoxels x Nparams x Nsubs
            sz = size(ifc_mintunks.img);
            ifc_cmro2 = copy(ifc_mintunks);
            ifc_cmro2.img = zeros(sz(1), sz(3));
            ifc_cmro2.fileprefix = strrep(ifc_cmro2.fileprefix, "mintunks", "cmro2");

            for idx_mg = 1:numel(mg)
                asub = this.parse_fileprefix(mg(idx_mg));
                ic_raichleks = this.raichleks( ...
                    sub=asub, ses="ses-*", input_func=opts.input_func, typeclass="nifti", stats="median");
                if isempty(ic_raichleks)
                    ic_raichleks = this.raichleks( ...
                        sub="sub-*", ses="ses-*", input_func=opts.input_func, typeclass="nifti", stats="median");
                end                
                ifc_cmro2.img(:,idx_mg) = ...
                    ifc_mintunks.img(:,1,idx_mg) .* ic_raichleks.imagingFormat.img(:,1) .* ...
                    60 * 44.64 * this.o2_content(asub);  % mmol / L / min
            end

            % apply opts.stats to img if img represents multiple mg;
            % do not apply opts.stats if there is only one mg
            if ~isemptytext(opts.stats) && ~isempty(ifc_cmro2.img) && numel(mg) > 1
                f = str2func(opts.stats);
                ifc_cmro2.img = f(ifc_cmro2.img, ndims(ifc_cmro2.img));  % apply opts.stats to Nses
                ifc_cmro2.fileprefix = ifc_cmro2.fileprefix + "-" + opts.stats;
            end

            ic = mlfourd.ImagingContext2(ifc_cmro2);
            obj = this.as(ic, opts.typeclass);            
        end

        function [obj,mg] = cmrglc(this, opts)
            %% mmol / L / min
            %  N.B.:  cmrglc := cmrglc * ; mmol glc / min / L := (10/180.156) mg glc / min / hg

            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            % for each available ic_mintunks:
            %     aufbau cmrglc from this.glc, ifc_martinv1, ic_huangks
            [ic_huangks,mg] = this.huangks(input_func=opts.input_func, typeclass="nifti", stats="");
            ifc_huangks = ic_huangks.imagingFormat;  % img ~ Nvoxels x Nparams x Nsubs
            sz = size(ifc_huangks.img);
            ifc_cmrglc = copy(ifc_huangks);
            ifc_cmrglc.img = zeros(sz(1), sz(3));
            ifc_cmrglc.fileprefix = strrep(ifc_cmrglc.fileprefix, "huangks", "cmrglc");

            for idx_mg = 1:numel(mg)
                asub = this.parse_fileprefix(mg(idx_mg));
                ic_martinv1 = this.martinv1( ...
                    sub=asub, ses="ses-*", input_func=opts.input_func, typeclass="nifti", stats="median");
                if isempty(ic_martinv1)
                    ic_martinv1 = this.martinv1( ...
                        sub="sub-*", ses="ses-*", input_func=opts.input_func, typeclass="nifti", stats="median");
                end 
                k1 = ifc_huangks.img(:,1,idx_mg);
                k2 = ifc_huangks.img(:,2,idx_mg);
                k3 = ifc_huangks.img(:,3,idx_mg);
                ifc_cmrglc.img(:,idx_mg) = ...
                    ic_martinv1.imagingFormat.img(:,1) .* k1 .* k3 ./ (k2 + k3);
                ifc_cmrglc.img(:,idx_mg) = ...
                    ifc_cmrglc.img(:,idx_mg) .* ...
                    ((1 / this.LC) * 60 * (10 / 180.156) * this.glc(asub));  % mmol / L / min
            end

            % apply opts.stats to img if img represents multiple mg;
            % do not apply opts.stats if there is only one mg
            if ~isemptytext(opts.stats) && ~isempty(ifc_cmrglc.img) && numel(mg) > 1
                f = str2func(opts.stats);
                ifc_cmrglc.img = f(ifc_cmrglc.img, ndims(ifc_cmrglc.img));  % apply opts.stats to Nses
                ifc_cmrglc.fileprefix = ifc_cmrglc.fileprefix + "-" + opts.stats;
            end

            ic = mlfourd.ImagingContext2(ifc_cmrglc);
            obj = this.as(ic, opts.typeclass);   
        end

        function obj = ogi(this, opts)
            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            [ic_cmro2,mg_cmro2] = this.cmro2(stats="", input_func=opts.input_func, typeclass="nifti");
            ifc_cmro2 = ic_cmro2.imagingFormat;
            [ic_cmrglc,mg_cmrglc] = this.cmrglc(stats="", input_func=opts.input_func, typeclass="nifti");
            ifc_cmrglc = ic_cmrglc.imagingFormat;
            img = [];

            % typically, only one FDG, but multiple O2 scans available per scan day
            for an_mg = asrow(mg_cmrglc)
                asub = this.parse_fileprefix(an_mg(1));
                select_cmro2 = contains(mg_cmro2, asub);  % as row
                select_cmrglc = contains(mg_cmrglc, asub);  % as row
                if sum(select_cmro2) >= 1 && sum(select_cmrglc) == 1
                    img_ = median(ifc_cmro2.img(:, select_cmro2), 2) ./ ifc_cmrglc.img(:, select_cmrglc);
                    img = [img, img_]; %#ok<AGROW>
                end
            end

            ifc_cmro2.img = img;
            ifc_cmro2.fileprefix = strrep(ifc_cmro2.fileprefix, "cmro2", "ogi");

            % apply opts.stats to img if img represents multiple mg;
            % do not apply opts.stats if there is only one mg
            if ~isemptytext(opts.stats) && ~isempty(ifc_cmro2.img) &&  size(ifc_cmro2.img, 2) > 1
                f = str2func(opts.stats);
                ifc_cmro2.img = f(ifc_cmro2.img, ndims(ifc_cmro2.img));  % apply opts.stats to Nses
                ifc_cmro2.fileprefix = ifc_cmro2.fileprefix + "-" + opts.stats;
            end

            ic = mlfourd.ImagingContext2(ifc_cmro2);
            obj = this.as(ic, opts.typeclass);
        end

        function obj = agi(this, opts)
            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            [ic_cmro2,mg_cmro2] = this.cmro2(stats="", input_func=opts.input_func, typeclass="nifti");
            ifc_cmro2 = ic_cmro2.imagingFormat;
            [ic_cmrglc,mg_cmrglc] = this.cmrglc(stats="", input_func=opts.input_func, typeclass="nifti");
            ifc_cmrglc = ic_cmrglc.imagingFormat;
            img = [];

            % typically, only one FDG, but multiple O2 scans available per scan day
            for an_mg = asrow(mg_cmrglc)
                asub = this.parse_fileprefix(an_mg(1));
                select_cmro2 = contains(mg_cmro2, asub);  % as row
                select_cmrglc = contains(mg_cmrglc, asub);  % as row
                if sum(select_cmro2) >= 1 && sum(select_cmrglc) == 1
                    img_ =  squeeze(ifc_cmrglc.img(:, select_cmrglc)) - median(ifc_cmro2.img(:, select_cmro2), 2) / 6;
                    img = [img, img_]; %#ok<AGROW>
                end
            end

            ifc_cmro2.img = img;
            ifc_cmro2.fileprefix = strrep(ifc_cmro2.fileprefix, "cmro2", "agi");

            % apply opts.stats to img if img represents multiple mg;
            % do not apply opts.stats if there is only one mg
            if ~isemptytext(opts.stats) && ~isempty(ifc_cmro2.img) &&  size(ifc_cmro2.img, 2) > 1
                f = str2func(opts.stats);
                ifc_cmro2.img = f(ifc_cmro2.img, ndims(ifc_cmro2.img));  % apply opts.stats to Nses
                ifc_cmro2.fileprefix = ifc_cmro2.fileprefix + "-" + opts.stats;
            end

            ic = mlfourd.ImagingContext2(ifc_cmro2);
            obj = this.as(ic, opts.typeclass);
        end

        %% modelled data objects,
        %  which access objects constructed from long computations of PycharmProjects/dynesty/idif2024

        function [obj,mg] = martinv1(this, opts)
            arguments
                this mlvg.Idif2024
                opts.sub {mustBeText} = "sub-*"  % "sub-*", "sub-108293", ...
                opts.ses {mustBeText} = "ses-*"  % "ses-*", "ses-20210421144815"
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = ""  % "median", "mean", "iqr", "std"
            end
            if contains(opts.input_func, "boxcar")
                opts.input_func = "idif";
            end
            if contains(opts.input_func, "twil")
                opts.input_func = "twilite";
            end

            %% recursion

            if strcmp(opts.sub, "sub-*") && strcmp(opts.ses, "ses-*")
                mg = mglob( ...
                    fullfile( ...
                        this.derivatives_path, opts.sub, opts.ses, "pet", ...
                        sprintf("*-%s_martinv1.nii.gz", opts.input_func)));
                assert(~isemptytext(mg))    
                if isemptytext(mg)
                    obj = [];
                    mg = [];
                    return
                end            

                % append img
                img = [];
                ic = [];
                fp = sprintf("sub-all_ses-all_trc-co_proc-schaefer-%s-martinv1", opts.input_func);
                for an_mg = asrow(mg)
                    [asub,ases] = this.parse_fileprefix(an_mg);
                    ic = this.martinv1(sub=asub, ses=ases, input_func=opts.input_func, typeclass="nifti");
                    img = [img, ic.imagingFormat.img]; %#ok<AGROW>
                end

                % apply opts.stats to img
                if ~isemptytext(opts.stats) && ~isempty(img)
                    f = str2func(opts.stats);
                    img = f(img, 2);
                    fp = fp + "-" + opts.stats;
                end

                % assemble final
                ifc = ic.imagingFormat;
                ifc.img = img;
                ifc.fqfp = fullfile(this.derivatives_path, fp);
                ic = mlfourd.ImagingContext2(ifc);
                obj = this.as(ic, opts.typeclass);
                return
            end

            %% base case

            [ic,mg] = this.schaefer_metric( ...
                sub=opts.sub, ses=opts.ses, trc="trc-co", metric=opts.input_func+"_martinv1");
            if isempty(mg)
                obj = [];
                return
            end
            if strcmp(opts.input_func, "idif")
                ic = ic ./ this.RECOVERY_COEFFICIENT;
            end
            s = split(ic.fileprefix, "proc-");
            ic.fileprefix = sprintf("%sproc-schaefer-%s-martinv1", s(1), opts.input_func);
            obj = this.as(ic, opts.typeclass);
        end

        function [obj,mg] = raichleks(this, opts)
            arguments
                this mlvg.Idif2024
                opts.sub {mustBeText} = "sub-*"  % "sub-*", "sub-108293", ...
                opts.ses {mustBeText} = "ses-*"  % "ses-*", "ses-20210421152358"
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = ""  % "median", "mean", "iqr", "std"
            end
            [obj,mg] = this.modelks( ...
                sub=opts.sub, ses=opts.ses, input_func=opts.input_func, typeclass=opts.typeclass, stats=opts.stats, ...
                model="Raichle1983", trc="trc-ho", fhandle=@this.raichleks, fname="raichleks");
        end

        function [obj,mg] = mintunks(this, opts)
            arguments
                this mlvg.Idif2024
                opts.sub {mustBeText} = "sub-*"  % "sub-*", "sub-108293", ...
                opts.ses {mustBeText} = "ses-*"  % "ses-*", "ses-20210421152358"
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = ""  % "median", "mean", "iqr", "std"
            end            
            [obj,mg] = this.modelks( ...
                sub=opts.sub, ses=opts.ses, input_func=opts.input_func, typeclass=opts.typeclass, stats=opts.stats, ...
                model="Mintun1984", trc="trc-oo", fhandle=@this.mintunks, fname="mintunks");
        end
        
        function [obj,mg] = huangks(this, opts)
            arguments
                this mlvg.Idif2024
                opts.sub {mustBeText} = "sub-*"  % "sub-*", "sub-108293", ...
                opts.ses {mustBeText} = "ses-*"  % "ses-*", "ses-20210421152358"
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = ""  % "median", "mean", "iqr", "std"
            end            
            [obj,mg] = this.modelks( ...
                sub=opts.sub, ses=opts.ses, input_func=opts.input_func, typeclass=opts.typeclass, stats=opts.stats, ...
                model="Huang1980", trc="trc-fdg", fhandle=@this.huangks, fname="huangks");
        end
        
        function [obj,mg] = ichiseks(this)
        end

        function [obj,mg] = logz(this)
        end
        
        function [obj,mg] = bayes_factors(this)
        end

        function [obj,mg] = information(this)
        end

        function [obj,mg] = residual(this)
        end

        %% pimary data objects,
        %  which access objects constructed from long computations on listmode, e.g., mlsiemens.BainMoCo2

        function [obj,mg] = brainmoco2_static(this, opts)
            arguments
                this mlvg.Idif2024
                opts.sub {mustBeTextScalar} = "sub-108293"
                opts.ses {mustBeTextScalar} = "ses-*"
                opts.trc {mustBeTextScalar} = "trc-ho"
            end

            mg = mglob( ...
                fullfile( ...
                    this.derivatives_path, opts.sub, opts.ses, "pet", ...
                     opts.sub + "_" + opts.ses + "_" + opts.trc + "_proc-delay0-BrainMoCo2-createNiftiStatic.nii.gz"));
            if isemptytext(mg)
                error("mlvg:FileNotFoundError", stackstr())
            end
            obj = mlfourd.ImagingContext2(mg(1));
        end
        
        function [obj,mg] = brainmoco2_maframes(this, opts)
            arguments
                this mlvg.Idif2024
                opts.sub {mustBeTextScalar} = "sub-108293"
                opts.ses {mustBeTextScalar} = "ses-*"
                opts.trc {mustBeTextScalar} = "trc-ho"
            end

            for tasi = 1:length(this.TIME_APPEND_SUFFIXES)
                tas = this.TIME_APPEND_SUFFIXES(tasi);
                mg = mglob( ...
                    fullfile( ...
                        this.sourcedata_path, opts.sub, opts.ses, "pet", ...
                        opts.sub + "_" + opts.ses + "_" + opts.trc + ...
                        "_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames" + tas + ".nii.gz"));
                if ~isemptytext(mg)
                    break
                end
            end
            if isemptytext(mg)
                error("mlvg:FileNotFoundError", stackstr())
            end
            obj = mlfourd.ImagingContext2(mg(1));
        end
        
        function [obj,mg] = brainmoco2_reshape_to_schaefer(this, opts)
            arguments
                this mlvg.Idif2024
                opts.sub {mustBeTextScalar} = "sub-108293"
                opts.ses {mustBeTextScalar} = "ses-*"
                opts.trc {mustBeTextScalar} = "trc-ho"
            end

            for tasi = 1:length(this.TIME_APPEND_SUFFIXES)
                tas = this.TIME_APPEND_SUFFIXES(tasi);
                mg = mglob( ...
                    fullfile( ...
                        this.derivatives_path, opts.sub, opts.ses, "pet", ...
                        opts.sub + "_" + opts.ses + "_" + opts.trc + ...
                        "_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames" + tas + ...
                        "-ParcSchaeffer-reshape-to-schaeffer-schaeffer.nii.gz"));
                if ~isemptytext(mg)
                    break
                end
            end
            if isemptytext(mg)
                error("mlvg:FileNotFoundError", stackstr())
            end
            obj = mlfourd.ImagingContext2(mg(1));      
        end
        
        function [obj,mg] = schaefer_metric(this, opts)
            arguments
                this mlvg.Idif2024
                opts.sub {mustBeTextScalar} = "sub-108293"
                opts.ses {mustBeTextScalar} = "ses-*"
                opts.trc {mustBeTextScalar} = "trc-co"
                opts.metric {mustBeTextScalar} = "idif_martinv1"
            end

            for tasi = 1:length(this.TIME_APPEND_SUFFIXES)
                tas = this.TIME_APPEND_SUFFIXES(tasi);
                mg = mglob( ...
                    fullfile( ...
                        this.derivatives_path, opts.sub, opts.ses, "pet", ...
                        opts.sub + "_" + opts.ses + "_" + opts.trc + ...
                        "_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames" + tas + ...
                        "*-schaeffer-" + opts.metric + ".nii.gz"));
                if ~isemptytext(mg)
                    break
                end
            end
            if isemptytext(mg)
                obj = [];
                mg = [];
                return
            end
            obj = mlfourd.ImagingContext2(mg(1));
        end

        %% construction

        function this = Idif2024(varargin)
        end
    end

    methods (Static)
        function o2 = calc_o2content(opts)
            %% CaO2 = Hb (gm/dl) x 1.34 ml O2/gm Hb x SaO2 + PaO2 x (.003 ml O2/mm Hg/dl) ~ mL/dL
            %  normal oxyhb ~ [11.7, 16.63], mean ~ 14.2
            %  normal SaO2 ~ [90, 95], mean ~ 92.5
            %  normal PaO2 ~ [83, 108], mean ~ 95.5

            arguments
                opts.oxyhb double = 14.2
                opts.SaO2 double = 0.925
                opts.PaO2 double = 95.5
            end
            if opts.SaO2 > 1
                opts.SaO2 = opts.SaO2/100;
            end

            o2 = 1.34 * opts.oxyhb * opts.SaO2 + 0.0031 * opts.PaO2;
        end
        
        function [sub,ses,trc] = parse_fileprefix(fp)
            re = regexp(fp, "(?<sub>sub-[A-Za-z0-9]+)_(?<ses>ses-[A-Za-z0-9]+)_trc\-(?<trc>[a-z]+)_", "names");
            sub = re.sub;
            ses = re.ses;
            trc = re.trc;
        end
    end

    %% PRIVATE

    properties (Access = private)
        glc_
        o2_content_
    end

    methods (Access = private)
        function obj1 = as(this, obj, typeclass)
            %  Args:  
            %      this mlvg.Idif2024
            %      obj {mustBeNonempty} : must be understood by mlfourd.ImagingContext2()
            %      typeclass {mustBeText} : 
            %          "nifti", "cifti", "fqfn", "double", "single"
            
            arguments
                this mlvg.Idif2024
                obj {mustBeNonempty}
                typeclass {mustBeText}
            end

            switch convertStringsToChars(typeclass)
                case {'nifti', 'nii.gz', 'nii'}
                    obj1 = obj;
                case {'cifti', 'dscalar.nii', 'dtseries.nii'}
                    obj1 = this.nifti_to_cifti(obj);
                case 'fqfn'
                    obj1 = obj.fqfn;
                case 'double'
                    obj1 = double(obj);
                case 'single'
                    obj1 = single(obj);
                otherwise
                    error("mlvg:ValueError", stackstr())
            end
        end
        
        function [obj,mg] = modelks(this, opts)
            %% does not apply opts.stats if there is only one mg

            arguments
                this mlvg.Idif2024
                opts.sub {mustBeText} = "sub-*"  % "sub-*", "sub-108293", ...
                opts.ses {mustBeText} = "ses-*"  % "ses-*", "ses-20210421152358"
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = ""  % "median", "mean", "iqr", "std"
                opts.model {mustBeTextScalar}
                opts.trc {mustBeTextScalar}
                opts.fhandle function_handle
                opts.fname {mustBeTextScalar}
            end
            if contains(opts.input_func, "boxcar")
                opts.input_func = "idif";
            end
            if contains(opts.input_func, "twil")
                opts.input_func = "twilite";
            end

            % metric
            if strcmp(opts.input_func, "idif")
                metric = opts.model + "Boxcar-main6-rc1p851000-qm";
            elseif strcmp(opts.input_func, "twilite")
                metric = opts.model + "Artery-main6-rc1p851000-qm";
            else
                error("mlvg:ValueError", stackstr())
            end

            %% recursion

            if strcmp(opts.sub, "sub-*") || strcmp(opts.ses, "ses-*")
                mg = mglob( ...
                    fullfile( ...
                        this.derivatives_path, opts.sub, opts.ses, "pet", ...
                        sprintf("*-%s.nii.gz", metric)));
                if isemptytext(mg)
                    obj = [];
                    mg = [];
                    return
                end

                % append img
                img = [];
                ic = [];
                fp = sprintf("sub-all_ses-all_%s_proc-schaefer-%s-%s", opts.trc, opts.input_func, opts.fname);
                for an_mg = asrow(mg)
                    [asub,ases] = this.parse_fileprefix(an_mg);
                    ic = opts.fhandle(sub=asub, ses=ases, input_func=opts.input_func, typeclass="nifti");
                    img = [img, ic.imagingFormat.img]; %#ok<AGROW>
                end
                img = reshape(img, [size(ic.imagingFormat.img), numel(mg)]);  % e.g., Nvoxels x Nparams x Nses

                % apply opts.stats to img if img represents multiple mg;
                % do not apply opts.stats if there is only one mg
                if ~isemptytext(opts.stats) && ~isempty(img) && numel(mg) > 1
                    f = str2func(opts.stats);
                    img = f(img, ndims(img));  % apply opts.stats to Nses
                    fp = fp + "-" + opts.stats;
                end

                % assemble final
                ifc = ic.imagingFormat;
                ifc.img = img;
                ifc.fqfp = fullfile(this.derivatives_path, fp);
                ic = mlfourd.ImagingContext2(ifc);
                obj = this.as(ic, opts.typeclass);
                return
            end

            %% base case
            [ic,mg] = this.schaefer_metric( ...
                sub=opts.sub, ses=opts.ses, trc=opts.trc, metric=metric);
            if isempty(mg)
                obj = [];
                return
            end
            s = split(ic.fileprefix, "proc-");
            ic.fileprefix = sprintf("%sproc-schaefer-%s-%s", s(1), opts.input_func, opts.fname);
            obj = this.as(ic, opts.typeclass);
        end

        function cii = nifti_to_cifti(this, nii)
            arguments
                this mlvg.Idif2024
                nii mlfourd.ImagingContext2
            end
            assert(contains(nii.fileprefix, "schaef"))

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
            img = nii.imagingFormat.img;
            if ndims(img) > 2 %#ok<ISMAT>
                sz = size(img);
                img = reshape(img, [sz(1), prod(sz(2:end))]);
            end
            Nt = size(img, 2);
            ld = load( ...                
                fullfile( ...
                    getenv("SINGULARITY_HOME"), ...
                    "CCIR_01211", "derivatives", "sub-108293", "ses-20210218", "Parcellations", ...
                    "indices_all.mat"));
            indices_all = ld.indices_all;  % no zero
            cdata1 = zeros(size(schaefer_indices, 1), Nt, "single");
            for i_ = asrow(unique_schaefer_indices)
                found = indices_all == i_;
                for t_ = 1:Nt
                    cdata1(i_ == schaefer_indices, t_) = img(found, t_);  % must be scalar
                end
            end

            % assemble final
            cii = schaefer;
            cii.cdata = cdata1;
            if size(cdata1, 2) > 1
                cii.diminfo{2} = cifti_diminfo_make_series(size(cdata1, 2), 0, 1, 'SECOND');
                cii.fqfn = convertStringsToChars(nii.fqfp + ".dtseries.nii");
            else
                cii.diminfo{2} = cifti_diminfo_make_scalars(1);
                cii.fqfn = convertStringsToChars(nii.fqfp + ".dscalar.nii");
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
