classdef Registration < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
	%% REGISTRATION  

	%  $Revision$
 	%  was created 10-Nov-2021 14:12:37 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
 	%% It was developed on Matlab 9.11.0.1769968 (R2021b) for MACI64.  Copyright 2021 John Joowon Lee.
 	
    methods (Static)
        function this = flirt_all_to_mni152(varargin)
            %  Returns:
            %      this: mlvg.Registration instance
        end
        function this = flirt_all_to_t1w(varargin)
            %  Args:
            %      projectPath (isfolder): e.g., fullfile(getenv("SINGULARITY_HOME"), "CCIR_01211"), 
            %                              needed by mlsiemens.BiographBids
            %      subjectFolder (isfolder): e.g., "sub-108293", 
            %                                needed by mlsiemens.BiographBids
            %  Returns:
            %      this: mlvg.Registration instance

            bids = mlvg.Ccir1211BidsBids(varargin{:});

            this = mlvg.Registration("bids", bids, varargin{:});            
            this.flirt_t2w_to_t1w();
            this.flirt_flair_to_t1w();
            this.flirt_tof_to_t1w();
            this.flirt_agtracers_to_t1w();
        end
        function this = flirt_pet_to_tof(varargin)
            %  Returns:
            %      this: mlvg.Registration instance
        end
    end

	properties
 		bids
        flirt
        reuse
        t1w_ic
        t1w_weight_ic
        tof_weight_ic
        wmparc_on_t1w_ic
 	end

	methods 
		  
 		function this = Registration(varargin)
            %  Params:
            %      bids (mlsiemens.BiographBids): 

            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, "bids", @(x) isa(x, "mlsiemens.BiographBids"))
            addParameter(ip, "reuse", true, @islogical)
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.reuse = ipr.reuse;
 			
            this.bids = copy(ipr.bids);
            this.t1w_ic = this.bids.t1w_ic;
            this.wmparc_on_t1w_ic = this.build_wmparc_on_t1w_ic();
            this.t1w_weight_ic = this.wmparc_on_t1w_ic.binarized();
            this.flirt = mlfsl.Flirt( ...
                "in", [], "inweight", [], ...
                "ref", this.t1w_ic.fqfn, "refweight", this.t1w_weight_ic.fqfn, ...
                "out", []);
        end

        function ic = flirt_t2w_to_t1w(this)
            flirt_ = copy(this.flirt);
            flirt_.in = this.bids.t2w_ic.fqfn;
            fn = strcat(this.bids.t2w_ic.fileprefix, "_on_T1w.nii.gz");
            flirt_.out = fullfile(this.bids.anatPath, fn);
            flirt_.flirt();
            ic = mlfourd.ImagingContext2(flirt_.out);
        end
        function ic = flirt_flair_to_t1w(this)
            flirt_ = copy(this.flirt);
            flirt_.in = this.bids.flair_ic.fqfn;
            fn = strcat(this.bids.flair_ic.fileprefix, "_on_T1w.nii.gz");
            flirt_.out = fullfile(this.bids.anatPath, fn);
            flirt_.flirt();
            ic = mlfourd.ImagingContext2(flirt_.out);
        end
        function ic = flirt_t1w_to_tof(this)
            
            % lazy tof_weight_ic
            if isempty(this.tof_weight_ic)
                this.tof_weight_ic = this.bids.tof_ic.thresh(20);
                this.tof_weight_ic = this.tof_weight_ic.binarized();
                this.tof_weight_ic.filepath = this.bids.anatPath;
                this.tof_weight_ic.save();
            end

%             out_xyztrans = fullfile(this.bids.anatPath, ...
%                 strcat(this.bids.t1w_ic.fileprefix, "_on_angio_xyztrans.nii.gz"));
%             flirt_ = mlfsl.Flirt( ...
%                 "in", this.bids.t1w_ic.fqfn, "inweight", "", ...
%                 "ref", this.bids.tof_ic.fqfn, "refweight", "", ...
%                 "schedule", fullfile(getenv("FSLDIR"), "etc", "flirtsch", "xyztrans.sch"), ...
%                 "cost", "normmi", ...
%                 "out", out_xyztrans);
%             flirt_.flirt();

            out = fullfile(this.bids.anatPath, ...
                strcat(this.bids.t1w_ic.fileprefix, "_on_angio.nii.gz"));
            flirt_ = mlfsl.Flirt( ...
                "in", this.bids.t1w_ic.fqfn, "inweight", this.t1w_weight_ic, ...
                "ref", this.bids.tof_ic.fqfn, "refweight", this.tof_weight_ic, ...
                "searchrx", 5, ...
                "cost", "normmi", ...
                "interp", "spline", ...
                "out", out);
            flirt_.flirt();

            ic = mlfourd.ImagingContext2(flirt_.out);
        end
        function ic = flirt_tof_to_t1w(this)

            % lazy t1w_on_tof
            out = fullfile(this.bids.anatPath, ...
                strcat(this.bids.t1w_ic.fileprefix, "_on_angio.nii.gz"));
            out_mat = strrep(out, '.nii.gz', '.mat');
            if ~isfile(out_mat) || ~this.reuse
                this.flirt_t1w_to_tof();
                assert(isfile(out_mat))
            end

            tof_on_t1w_ic = mlfourd.ImagingContext2( ...
                fullfile(this.bids.anatPath, ...
                strcat(this.bids.tof_ic.fileprefix, "_on_T1w.nii.gz")));
            flirt_ = mlfsl.Flirt( ...
                "in", this.bids.tof_ic.fqfn, "inweight", this.tof_weight_ic.fqfn, ...
                "ref", this.bids.t1w_ic.fqfn, "refweight", this.t1w_weight_ic.fqfn, ...
                "out", tof_on_t1w_ic.fqfn, ...
                "omat", out_mat);
            flirt_.invertXfm(); % omat -> init
            flirt_.applyXfm(); % init
            
            ic = mlfourd.ImagingContext2(flirt_.out);
        end
        function ics = flirt_agtracers_to_t1w(this)
        end

        function flirt_T1w_to_highest_res(this)
            %% flirt to 0.26 x 0.26 x 0.5 

            flirt = '/usr/local/fsl/bin/flirt';
            in = '/Users/jjlee/Singularity/CCIR_01211/sourcedata/sub-108293/anat/sub-108293_20210218081030_T1w.nii.gz';
            ref = '/Users/jjlee/Singularity/CCIR_01211/derivatives/sub-108293/anat/sub-108293_20210218081030_angio_embedded.nii.gz';
            out = '/Users/jjlee/Singularity/CCIR_01211/derivatives/sub-108293/anat/sub-108293_20210218081030_T1w_embedded.nii.gz';
            omat = '/Users/jjlee/Singularity/CCIR_01211/derivatives/sub-108293/anat/sub-108293_20210218081030_T1w_embedded.mat';
            refweight = '/Users/jjlee/Singularity/CCIR_01211/derivatives/sub-108293/anat/sub-108293_20210218081030_angio_embedded_binarized.nii.gz';
            inweight = '/Users/jjlee/Singularity/CCIR_01211/derivatives/sub-108293/anat/sub-108293_20210218081030_T1w_brain_mask.nii.gz';
            bins = 512;
            searchr = 5;
            cmd = sprintf( ...
                '%s -in %s -ref %s -out %s -omat %s -bins %i -cost normmi -searchrx -%i %i -searchry -%i %i -searchrz -%i %i -dof 6 -refweight %s -inweight %s -interp spline', ...
                flirt, in, ref, out, omat, bins, searchr, searchr, searchr, searchr, searchr, searchr, refweight, inweight);
            system(cmd)
        end
        function flirt_CT_to_highest_res(this)
            %% flirt to 0.26 x 0.26 x 0.5 

            flirt = '/usr/local/fsl/bin/flirt';
            in = '/Users/jjlee/Singularity/CCIR_01211/sourcedata/sub-108293/pet/sub-108293_20210421134537_CT_Brain_3.0_J30f_1_3.nii.gz';
            ref = '/Users/jjlee/Singularity/CCIR_01211/derivatives/sub-108293/anat/sub-108293_20210218081030_T1w_embedded.nii.gz';
            out = '/Users/jjlee/Singularity/CCIR_01211/derivatives/sub-108293/anat/sub-108293_20210421134537_CT_Brain_3.0_J30f_1_3_embedded.nii.gz';
            omat = '/Users/jjlee/Singularity/CCIR_01211/derivatives/sub-108293/anat/sub-108293_20210421134537_CT_Brain_3.0_J30f_1_3_embedded.mat';
            refweight = '/Users/jjlee/Singularity/CCIR_01211/derivatives/sub-108293/anat/sub-108293_20210218081030_T1w_embedded_binarized.nii.gz';
            inweight = '/Users/jjlee/Singularity/CCIR_01211/derivatives/sub-108293/anat/sub-108293_20210421134537_CT_Brain_3.0_J30f_1_3_thr0_binarized.nii.gz';
            bins = 512;
            searchr = 5;
            cmd = sprintf( ...
                '%s -in %s -ref %s -out %s -omat %s -bins %i -cost normmi -searchrx -%i %i -searchry -%i %i -searchrz -%i %i -dof 6 -refweight %s -inweight %s -interp spline', ...
                flirt, in, ref, out, omat, bins, searchr, searchr, searchr, searchr, searchr, searchr, refweight, inweight);
            system(cmd)
        end
        function flirt_xfm_to_highest_res(this)
            %% flirt to 0.26 x 0.26 x 0.5 

            xfm = 'sub-108293_20210218081030_T1w_embedded.mat';
        end
    end 

    %% PROTECTED

    methods (Access = protected)
        function ic = build_wmparc_on_t1w_ic(this)
            fn = strcat("wmparc_on_", this.bids.t1w_ic.filename);
            out = fullfile(this.bids.mriPath, fn);
            if isfile(out) && this.reuse
                ic = mlfourd.ImagingContext2(out);
                return
            end

            flirt_ = mlfsl.Flirt( ...
                "in", this.bids.T1_ic.fqfn, ...
                "ref", this.bids.t1w_ic.fqfn);
            flirt_.flirt();
            flirt_.in = this.bids.wmparc_ic.fqfn;
            flirt_.out = out;
            flirt_.interp = 'nearestneighbour';
            flirt_.applyXfm()
            ic = mlfourd.ImagingContext2(flirt_.out);
        end
        function that = copyElement(this)
            that = copyElement@matlab.mixin.Copyable(this);
            that.bids = copy(this.bids);
            that.flirt = copy(this.flirt);
            that.t1w_ic = copy(this.t1w_ic);
            that.t1w_weight_ic = copy(this.t1w_weight_ic);
            that.wmparc_on_t1w_ic = copy(this.wmparc_on_t1w_ic);
        end
    end    

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

