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

    properties
        fig_position = [80 80 1618 1000]; % coordinates for figures ~ [x0, y0, Dx, Dy]
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

        %% visualizations

        function [f,h] = rm_raincloud(this, opts)
            arguments
                this mlvg.Idif2024 
                opts.metric {mustBeTextScalar} = "CMRglc"
                opts.fname {mustBeTextScalar} = ""
                opts.axis_label {mustBeTextScalar} = "CMRglc (mmol L^{-1} min^{-1})"
                opts.idx_metric double = 1
                opts.rescale double = 1 
            end
            if isemptytext(opts.fname)
                opts.fname = lower(opts.metric);
            end

            %% colormapping

            [cb] = cbrewer2('Spectral', 12, 'pchip');
            cl(1, :) = cb(11, :);  % blue is IDIF
            cl(2, :) = cb(2, :);  % red is Twilite

            pwd0 = pushd(this.derivatives_path);

            %% aufbau data

            f = str2func(opts.fname);
            [nii_idif,mg_idif] = f(this, input_func="idif", stats="");
            [nii_twil,mg_twil] = f(this, input_func="twil", stats="");
            if ~isempty(getenv("DEBUG"))
                imagesc(nii_idif);
                imagesc(nii_twil);
            end

            % twilite data is more limited than idif data; use this.match_globbed()
            ifc_idif = nii_idif.imagingFormat;
            ifc_twil = nii_twil.imagingFormat;
            [ifc_idif, ifc_twil] = this.match_globbed(ifc_idif, mg_idif, ifc_twil, mg_twil);

            % select index of metric when necessary
            if ndims(ifc_idif) > 2 %#ok<ISMAT>
                ifc_idif.img = squeeze(ifc_idif.img(:,opts.idx_metric,:));
            end
            if ndims(ifc_twil) > 2 %#ok<ISMAT>
                ifc_twil.img = squeeze(ifc_twil.img(:,opts.idx_metric,:));
            end

            assert(all(size(ifc_idif.img) == size(ifc_twil.img)))
            N_sub = size(ifc_idif.img, ndims(ifc_idif.img));
            data = cell(N_sub, 2);  % N_sub x N. input func. methods
            for idx_sub = 1:N_sub
                data{idx_sub, 1} = opts.rescale * double(ifc_idif.img(:, idx_sub));
                data{idx_sub, 2} = opts.rescale * double(ifc_twil.img(:, idx_sub));
            end

            %% make & save figure

            f_metric  = figure(Position=this.fig_position);
            h = rm_raincloud(data, cl);
            %set(gca, 'YLim', [-0.3 1.6]);
            xlabel(opts.axis_label);
            ylabel("Participants");
            fontsize(scale=2.5)
            saveFigure2(f_metric, fullfile(pwd, lower(opts.metric) + "_raincloud"));

            popd(pwd0);
        end

        %% physiological data objects

        function [obj,mg] = cbv(this, opts)
            %% mL/mL median over all available data

            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            [ic,mg] = this.martinv1(input_func=opts.input_func, typeclass="nifti", stats=opts.stats);
            ifc = ic.imagingFormat;
            ifc.img = ifc.img ./ 0.85;  % R := ratio of small-vessel to large-vessel hematocrit, described by Mintun 1984
            ifc.fileprefix = strrep(ifc.fileprefix, "martinv1", "cbv");
            ic = mlfourd.ImagingContext2(ifc);
            obj = this.as(ic, opts.typeclass);
        end

        function [obj,mg] = cbf(this, opts)
            %% mL/mL/min median over all available data

            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            [ic,mg] = this.raichleks( ...
                idx_metric=1, ...
                input_func=opts.input_func, typeclass="nifti", stats=opts.stats);
            ic = ic * 60;
            ic.fileprefix = strrep(ic.fileprefix, "raichleks", "cbf");
            obj = this.as(ic, opts.typeclass);
        end

        function [obj,mg] = lambda(this, opts)
            %% partition coefficient for water, mL/mL

            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            [ic,mg] = this.raichleks( ...
                idx_metric=2, ...
                input_func=opts.input_func, typeclass="nifti", stats=opts.stats);
            ic.fileprefix = strrep(ic.fileprefix, "raichleks", "lambda");
            obj = this.as(ic, opts.typeclass);
        end

        function [obj,mg] = ps(this, opts)
            %% permeability-surface area product, mL/mL/min

            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            [ic,mg] = this.raichleks( ...
                idx_metric=3, ...
                input_func=opts.input_func, typeclass="nifti", stats=opts.stats);
            ic = ic * 60;
            ic.fileprefix = strrep(ic.fileprefix, "raichleks", "ps");
            obj = this.as(ic, opts.typeclass);
        end

        function [obj,mg] = oef(this, opts)
            %% mL/mL

            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            [ic,mg] = this.mintunks( ...
                idx_metric=1, ...
                input_func=opts.input_func, typeclass="nifti", stats=opts.stats);
            ic.fileprefix = strrep(ic.fileprefix, "mintunks", "oef");
            obj = this.as(ic, opts.typeclass);
        end

        function [obj,mg] = vcapillary(this, opts)
            %% r"$v_{post} + 0.5 v_{cap}$", mL/mL

            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            [ic,mg] = this.mintunks( ...
                idx_metric=2, ...
                input_func=opts.input_func, typeclass="nifti", stats=opts.stats);
            ic.fileprefix = strrep(ic.fileprefix, "mintunks", "vcapillary");
            obj = this.as(ic, opts.typeclass);
        end

        function [obj,mg] = fwatermetab(this, opts)
            %% r"frac. water of metab.", mL/mL

            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            [ic,mg] = this.mintunks( ...
                idx_metric=3, ...
                input_func=opts.input_func, typeclass="nifti", stats=opts.stats);
            ic.fileprefix = strrep(ic.fileprefix, "mintunks", "fwatermetab");
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

        function [obj,mg] = upper_k1(this, varargin)
            [obj,mg] = this.K1(varargin{:});
        end

        function [obj,mg] = K1(this, opts)
            %% mL / cm^3 / min

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
            ifc_K1 = copy(ifc_huangks);
            ifc_K1.img = zeros(sz(1), sz(3));
            ifc_K1.fileprefix = strrep(ifc_K1.fileprefix, "huangks", "K1");

            for idx_mg = 1:numel(mg)
                asub = this.parse_fileprefix(mg(idx_mg));
                ic_martinv1 = this.martinv1( ...
                    sub=asub, ses="ses-*", input_func=opts.input_func, typeclass="nifti", stats="median");
                if isempty(ic_martinv1)
                    ic_martinv1 = this.martinv1( ...
                        sub="sub-*", ses="ses-*", input_func=opts.input_func, typeclass="nifti", stats="median");
                end 
                k1 = ifc_huangks.img(:,1,idx_mg);
                ifc_K1.img(:,idx_mg) = ...
                    ic_martinv1.imagingFormat.img(:,1) .* k1 .* 60;  % mL / cm^3 / min
            end

            % apply opts.stats to img if img represents multiple mg;
            % do not apply opts.stats if there is only one mg
            if ~isemptytext(opts.stats) && ~isempty(ifc_K1.img) && numel(mg) > 1
                f = str2func(opts.stats);
                ifc_K1.img = f(ifc_K1.img, ndims(ifc_K1.img));  % apply opts.stats to Nses
                ifc_K1.fileprefix = ifc_K1.fileprefix + "-" + opts.stats;
            end

            ic = mlfourd.ImagingContext2(ifc_K1);
            obj = this.as(ic, opts.typeclass);   
        end

        function [obj,mg] = k2(this, opts)
            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            [obj,mg] = this.kss( ...
                idx_metric=2, ...
                input_func=opts.input_func, typeclass=opts.typeclass, stats=opts.stats);
        end

        function [obj,mg] = k3(this, opts)
            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            [obj,mg] = this.kss( ...
                idx_metric=3, ...
                input_func=opts.input_func, typeclass=opts.typeclass, stats=opts.stats);
        end

        function [obj,mg] = k4(this, opts)
            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
            end

            [obj,mg] = this.kss(  ...
                idx_metric=4, ...
                input_func=opts.input_func, typeclass=opts.typeclass, stats=opts.stats);
        end

        function [obj,mg] = kss(this, opts)
            %% k2 .. k4, in mL / cm^3

            arguments
                this mlvg.Idif2024
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = "median"  % "median", "mean", "iqr", "std"
                opts.idx_metric double = 2
            end

            [ic_huangks,mg] = this.huangks(input_func=opts.input_func, typeclass="nifti", stats="");
            ifc_huangks = ic_huangks.imagingFormat;  % img ~ Nvoxels x Nparams x Nsubs
            sz = size(ifc_huangks.img);
            ifc_kss = copy(ifc_huangks);
            ifc_kss.img = zeros(sz(1), sz(3));
            ifc_kss.fileprefix = strrep(ifc_kss.fileprefix, "huangks", "k" + opts.idx_metric);

            ifc_kss.img = ifc_huangks.img(:,opts.idx_metric,:) .* 60;  % mL / cm^3

            % apply opts.stats to img if img represents multiple mg;
            % do not apply opts.stats if there is only one mg
            if ~isemptytext(opts.stats) && ~isempty(ifc_kss.img) && numel(mg) > 1
                f = str2func(opts.stats);
                ifc_kss.img = f(ifc_kss.img, ndims(ifc_kss.img));  % apply opts.stats to Nses
                ifc_kss.fileprefix = ifc_kss.fileprefix + "-" + opts.stats;
            end

            ic = mlfourd.ImagingContext2(ifc_kss);
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

        function [obj,mg_cmrglc] = ogi(this, opts)
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

            % typically, only one FDG, but multiple O2 scans available per scan day
            img = [];
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

        function [obj,mg_cmrglc] = agi(this, opts)
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

            % typically, only one FDG, but multiple O2 scans available per scan day
            img = [];
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
                opts.model {mustBeTextScalar} = "Raichle1983"
                opts.trc {mustBeTextScalar} = "trc-ho"
                opts.idx_metric double = []
            end
            [obj,mg] = this.modelks( ...
                sub=opts.sub, ses=opts.ses, input_func=opts.input_func, typeclass=opts.typeclass, stats=opts.stats, ...
                model=opts.model, trc=opts.trc, fhandle=@this.raichleks, fname="raichleks");
            ifc = obj.imagingFormat;
            if ~isempty(opts.idx_metric)
                if isemptytext(opts.stats)
                    ifc.img = ifc.img(:, opts.idx_metric, :);
                else
                    ifc.img = ifc.img(:, opts.idx_metric);
                end
            end
            obj = mlfourd.ImagingContext2(ifc);
        end

        function [obj,mg] = mintunks(this, opts)
            arguments
                this mlvg.Idif2024
                opts.sub {mustBeText} = "sub-*"  % "sub-*", "sub-108293", ...
                opts.ses {mustBeText} = "ses-*"  % "ses-*", "ses-20210421152358"
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = ""  % "median", "mean", "iqr", "std"
                opts.model {mustBeTextScalar} = "Mintun1984"
                opts.trc {mustBeTextScalar} = "trc-oo"
                opts.idx_metric double = []
            end            
            [obj,mg] = this.modelks( ...
                sub=opts.sub, ses=opts.ses, input_func=opts.input_func, typeclass=opts.typeclass, stats=opts.stats, ...
                model=opts.model, trc=opts.trc, fhandle=@this.mintunks, fname="mintunks");
            ifc = obj.imagingFormat;
            if ~isempty(opts.idx_metric)
                if isemptytext(opts.stats)
                    ifc.img = ifc.img(:, opts.idx_metric, :);
                else
                    ifc.img = ifc.img(:, opts.idx_metric);
                end
            end
            obj = mlfourd.ImagingContext2(ifc);
        end
        
        function [obj,mg] = huangks(this, opts)
            arguments
                this mlvg.Idif2024
                opts.sub {mustBeText} = "sub-*"  % "sub-*", "sub-108293", ...
                opts.ses {mustBeText} = "ses-*"  % "ses-*", "ses-20210421152358"
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = ""  % "median", "mean", "iqr", "std"
                opts.model {mustBeTextScalar} = "Huang1980"
                opts.trc {mustBeTextScalar} = "trc-fdg"
                opts.idx_metric double = []
            end            
            [obj,mg] = this.modelks( ...
                sub=opts.sub, ses=opts.ses, input_func=opts.input_func, typeclass=opts.typeclass, stats=opts.stats, ...
                model=opts.model, trc=opts.trc, fhandle=@this.huangks, fname="huangks");
            ifc = obj.imagingFormat;
            if ~isempty(opts.idx_metric)
                if isemptytext(opts.stats)
                    ifc.img = ifc.img(:, opts.idx_metric, :);
                else
                    ifc.img = ifc.img(:, opts.idx_metric);
                end
            end
            obj = mlfourd.ImagingContext2(ifc);
        end
        
        function [obj,mg] = ichiseks(this)
        end
        
        function [obj,mg] = bayes_factors(this, opts)
            arguments
                this mlvg.Idif2024
                opts.sub {mustBeText} = "sub-*"  % "sub-*", "sub-108293", ...
                opts.ses {mustBeText} = "ses-*"  % "ses-*", "ses-20210421152358"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = ""  % "median", "mean", "iqr", "std"
                opts.model {mustBeTextScalar} = "Raichle1983"  
                opts.trc {mustBeTextScalar} = "trc-ho"  
            end            

            [ic_idif,mg_idif] = this.logz( ...
                sub=opts.sub, ses=opts.ses, input_func="idif", typeclass="nifti", stats="", ...
                model=opts.model, trc=opts.trc);  % retain sessions
            ifc_idif = ic_idif.imagingFormat;
            [ic_twil,mg_twil] = this.logz( ...
                sub=opts.sub, ses=opts.ses, input_func="twil", typeclass="nifti", stats="", ...
                model=opts.model, trc=opts.trc);  % retain sessions
            ifc_twil = ic_twil.imagingFormat;

            [ifc_idif, ifc_twil] = this.match_globbed(ifc_idif, mg_idif, ifc_twil, mg_twil);

            img = ifc_idif.img - ifc_twil.img;
            img = squeeze(img);  % awkwardness arises from this.modelks expecting vector params, but logz is scalar
            ifc_idif.img = img;
            ifc_idif.fileprefix = strrep(ifc_idif.fileprefix, "logz", "bfactors");

            % apply opts.stats to img if img represents multiple mg;
            % do not apply opts.stats if there is only one mg
            if ~isemptytext(opts.stats) && ~isempty(ifc_idif.img) && size(ifc_idif.img, 2) > 1
                f = str2func(opts.stats);
                ifc_idif.img = f(ifc_idif.img, ndims(ifc_idif.img));  % apply opts.stats to Nses
                ifc_idif.fileprefix = ifc_idif.fileprefix + "-" + opts.stats;
            end

            ic = mlfourd.ImagingContext2(ifc_idif);
            obj = this.as(ic, opts.typeclass);
        end

        function [obj,mg] = logz_ho(this, opts)
            arguments
                this mlvg.Idif2024
                opts.sub {mustBeText} = "sub-*"  % "sub-*", "sub-108293", ...
                opts.ses {mustBeText} = "ses-*"  % "ses-*", "ses-20210421152358"
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = ""  % "median", "mean", "iqr", "std"
            end

            [obj,mg] = this.logz( ...
                sub=opts.sub, ses=opts.ses, input_func=opts.input_func, typeclass=opts.typeclass, stats=opts.stats, ...
                model="Raichle1983", trc="trc-ho");
        end

        function [obj,mg] = logz_oo(this, opts)
            arguments
                this mlvg.Idif2024
                opts.sub {mustBeText} = "sub-*"  % "sub-*", "sub-108293", ...
                opts.ses {mustBeText} = "ses-*"  % "ses-*", "ses-20210421152358"
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = ""  % "median", "mean", "iqr", "std"
            end  

            [obj,mg] = this.logz( ...
                sub=opts.sub, ses=opts.ses, input_func=opts.input_func, typeclass=opts.typeclass, stats=opts.stats, ...
                model="Mintun1984", trc="trc-oo");
        end

        function [obj,mg] = logz_fdg(this, opts)
            arguments
                this mlvg.Idif2024
                opts.sub {mustBeText} = "sub-*"  % "sub-*", "sub-108293", ...
                opts.ses {mustBeText} = "ses-*"  % "ses-*", "ses-20210421152358"
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = ""  % "median", "mean", "iqr", "std"
            end  

            [obj,mg] = this.logz( ...
                sub=opts.sub, ses=opts.ses, input_func=opts.input_func, typeclass=opts.typeclass, stats=opts.stats, ...
                model="Huang1980", trc="trc-fdg");
        end

        function [obj,mg] = logz(this, opts)
            arguments
                this mlvg.Idif2024
                opts.sub {mustBeText} = "sub-*"  % "sub-*", "sub-108293", ...
                opts.ses {mustBeText} = "ses-*"  % "ses-*", "ses-20210421152358"
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = ""  % "median", "mean", "iqr", "std"
                opts.model {mustBeTextScalar} = "Raichle1983"  
                opts.trc {mustBeTextScalar} = "trc-ho"  
            end            
            [obj,mg] = this.modelks( ...
                sub=opts.sub, ses=opts.ses, input_func=opts.input_func, typeclass=opts.typeclass, stats=opts.stats, ...
                model=opts.model, trc=opts.trc, fhandle=@this.logz, fname="logz", ...
                metric="logz");
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
        function obj1 = as(obj, typeclass)
            %  Args:  
            %      this mlvg.Idif2024
            %      obj {mustBeNonempty} : must be understood by mlfourd.ImagingContext2()
            %      typeclass {mustBeText} : 
            %          "nifti", "cifti", "fqfn", "double", "single"
            
            arguments
                obj {mustBeNonempty}
                typeclass {mustBeText}
            end

            switch convertStringsToChars(typeclass)
                case {'nifti', 'nii.gz', 'nii'}
                    obj1 = obj;
                case {'cifti', 'dscalar.nii', 'dtseries.nii'}
                    obj1 = mlvg.Idif2024.nifti_to_cifti(obj);
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

        function cii = nifti_to_cifti(nii)
            arguments
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
        function [ifc_idif, ifc_twil] = match_globbed(this, ifc_idif, mg_idif, ifc_twil, mg_twil)
            %% There may be more boxcar results than twilite results:  N_idif >= N_twil.
            %  Twilite results are injective into boxcar results.

            mg_twil1 = strrep(mg_twil, "Artery", "Boxcar");
            mg_twil1 = strrep(mg_twil1, "twilite", "idif");
            matching = mg_idif' == mg_twil1;  % N_idif x N_twil
            matching = logical(sum(matching, 2))';  % 1 x N_idif

            if ndims(ifc_idif) == 2 %#ok<ISMAT>
                ifc_idif.img = ifc_idif.img(:,matching);
            end
            if ndims(ifc_idif) == 3
                ifc_idif.img = ifc_idif.img(:,:,matching);
            end
            if ndims(ifc_idif) == 4
                ifc_idif.img = ifc_idif.img(:,:,:,matching);
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
                opts.fname {mustBeTextScalar}  % corresponding to opts.fhandle, e.g., "raichleks", "logz", ...
                opts.metric {mustBeTextScalar} = "qm"  % "qm", "logz", "information", "residual"
            end
            if contains(opts.input_func, "boxcar")
                opts.input_func = "idif";
            end
            if contains(opts.input_func, "twil")
                opts.input_func = "twilite";
            end

            % metric
            if strcmp(opts.input_func, "idif")
                metric = opts.model + "Boxcar-main6-rc1p851000-" + opts.metric;
            elseif strcmp(opts.input_func, "twilite")
                metric = opts.model + "Artery-main6-rc1p851000-" + opts.metric;
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
                if strcmp(opts.sub, "sub-*")
                    sub_ = "sub-all";
                else
                    sub_ = opts.sub;
                end
                if strcmp(opts.ses, "ses-*")
                    ses_ = "ses-all";
                else
                    ses_ = opts.ses;
                end
                fp = sprintf("%s_%s_%s_proc-schaefer-%s-%s", sub_, ses_, opts.trc, opts.input_func, opts.fname);
                for an_mg = asrow(mg)
                    [asub,ases] = this.parse_fileprefix(an_mg);
                    ic = opts.fhandle( ...
                        sub=asub, ses=ases, input_func=opts.input_func, typeclass="nifti", ...
                        model=opts.model, trc=opts.trc);
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
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
