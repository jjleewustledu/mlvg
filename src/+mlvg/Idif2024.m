classdef Idif2024 < handle & mlsystem.IHandle
    %% Provides simpler, high-level, re-usable application programming interfaces for software supporting
    %  Lee, et al., to be submitted to Physics in Biology & Medicine, 2024.
    %
    %  See also:  mlvg.Lee2024, PycharmProjects/dynesty/idif2024, Singularity/CCIR_01222.
    %  
    %  Created 06-Mar-2024 23:13:15 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 23.2.0.2515942 (R2023b) Update 7 for MACA64.  Copyright 2024 John J. Lee.
    

    properties (Constant)
        RECOVERY_COEFFICIENT = 1.8509  % determined from CO for Biograph Vision 600
        TIME_APPEND_SUFFIXES = ["_timeAppend-165", "_timeAppend-80", "_timeAppend-4", ""]
    end

    properties (Dependent)
        derivatives_path
        o2_content  % containers.map:  sub-id =: o2 content ~ [0, 1] in mL/mL
        sourcedata_path
    end

    methods  %% GET
        function g = get.derivatives_path(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "CCIR_01211", "derivatives");
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

        function obj = cbv(this)
            %% mL/mL
        end

        function obj = cbf(this)
            %% mL/mL/min
        end

        function obj = lambda(this)
            %% partition coefficient for water, mL/mL
        end

        function obj = ps(this)
            %% permeability-surface area product, mL/mL/min
        end

        function obj = oef(this)
        end

        function obj = v_capillary(this)
            %% r"$v_{post} + 0.5 v_{cap}$", mL/mL
        end

        function obj = f_water_metab(this)
            %% r"frac. water of metab."
        end

        function obj = cmro2(this)
            %% N.B.:  cmro2 := cmro2 * 44.64; \mu mol O_2 := mL O_2

        end

        function obj = cmrglc(this)
        end

        function obj = ogi(this)
        end

        function obj = agi(this)
        end

        %% modelled data objects,
        %  which access objects constructed from long computations of PycharmProjects/dynesty/idif2024

        function obj = martinv1(this, opts)
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

            ic = this.schaefer_metric( ...
                sub=opts.sub, ses=opts.ses, trc="trc-co", metric=opts.input_func+"_martinv1");
            if strcmp(opts.input_func, "idif")
                ic = ic ./ this.RECOVERY_COEFFICIENT;
            end
            s = split(ic.fileprefix, "proc-");
            ic.fileprefix = sprintf("%sproc-schaefer-%s-martinv1", s(1), opts.input_func);
            obj = this.as(ic, opts.typeclass);
        end

        function obj = raichleks(this, opts)
            arguments
                this mlvg.Idif2024
                opts.sub {mustBeText} = "sub-*"  % "sub-*", "sub-108293", ...
                opts.ses {mustBeText} = "ses-*"  % "ses-*", "ses-20210421152358"
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

            % metric
            if strcmp(opts.input_func, "idif")
                metric = "Raichle1983Boxcar-main6-rc1p851000-qm";
            elseif strcmp(opts.input_func, "twilite")
                metric = "Raichle1983Artery-main6-rc1p851000-qm";
            else
                error("mlvg:ValueError", stackstr())
            end

            %% recursion

            if strcmp(opts.sub, "sub-*") && strcmp(opts.ses, "ses-*")
                mg = mglob( ...
                    fullfile( ...
                        this.derivatives_path, opts.sub, opts.ses, "pet", ...
                        sprintf("*-%s.nii.gz", metric)));
                assert(~isemptytext(mg))                

                % append img
                img = [];
                ic = [];
                fp = sprintf("sub-all_ses-all_trc-ho_proc-schaefer-%s-raichleks", opts.input_func);
                for an_mg = asrow(mg)
                    [asub,ases] = this.parse_fileprefix(an_mg);
                    ic = this.raichleks(sub=asub, ses=ases, input_func=opts.input_func, typeclass="nifti");
                    img = [img, ic.imagingFormat.img]; %#ok<AGROW>
                end
                img = reshape(img, [size(ic.imagingFormat.img), numel(mg)]);

                % apply opts.stats to img
                if ~isemptytext(opts.stats) && ~isempty(img)
                    f = str2func(opts.stats);
                    img = f(img, ndims(img));
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
            ic = this.schaefer_metric( ...
                sub=opts.sub, ses=opts.ses, trc="trc-ho", metric=metric);
            s = split(ic.fileprefix, "proc-");
            ic.fileprefix = sprintf("%sproc-schaefer-%s-raichleks", s(1), opts.input_func);
            obj = this.as(ic, opts.typeclass);
        end

        function obj = mintunks(this, opts)
            arguments
                this mlvg.Idif2024
                opts.sub {mustBeText} = "sub-*"  % "sub-*", "sub-108293", ...
                opts.ses {mustBeText} = "ses-*"  % "ses-*", "ses-20210421152358"
                opts.input_func {mustBeTextScalar} = "idif"  % "idif", "boxcar", "aif", "twil", "twilite"
                opts.typeclass {mustBeTextScalar} = "nifti"  % "nifti", "cifti"
                opts.stats {mustBeTextScalar} = ""  % "median", "mean", "iqr", "std"
            end            
            copts = namedargs2cell(opts);
            obj = this.modelks( ...
                sub=opts.sub, ses=opts.ses, input_func=opts.input_func, typeclass=opts.typeclass, stats=opts.stats, ...
                model="Mintun1984", trc="trc-oo", fhandle=@this.mintunks, fname="mintunks");
        end
        
        function obj = huangks(this)
        end
        
        function obj = ichiseks(this)
        end

        function obj = logz(this)
        end
        
        function obj = bayes_factors(this)
        end

        function obj = information(this)
        end

        function obj = residual(this)
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
                error("mlvg:FileNotFoundError", stackstr())
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
        
        function obj = modelks(this, opts)
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

            if strcmp(opts.sub, "sub-*") && strcmp(opts.ses, "ses-*")
                mg = mglob( ...
                    fullfile( ...
                        this.derivatives_path, opts.sub, opts.ses, "pet", ...
                        sprintf("*-%s.nii.gz", metric)));
                assert(~isemptytext(mg))                

                % append img
                img = [];
                ic = [];
                fp = sprintf("sub-all_ses-all_%s_proc-schaefer-%s-mintunks", opts.trc, opts.input_func);
                for an_mg = asrow(mg)
                    [asub,ases] = this.parse_fileprefix(an_mg);
                    ic = opts.fhandle(sub=asub, ses=ases, input_func=opts.input_func, typeclass="nifti");
                    img = [img, ic.imagingFormat.img]; %#ok<AGROW>
                end
                img = reshape(img, [size(ic.imagingFormat.img), numel(mg)]);

                % apply opts.stats to img
                if ~isemptytext(opts.stats) && ~isempty(img)
                    f = str2func(opts.stats);
                    img = f(img, ndims(img));
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
            ic = this.schaefer_metric( ...
                sub=opts.sub, ses=opts.ses, trc=opts.trc, metric=metric);
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
