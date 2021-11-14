classdef Registration < handle
	%% REGISTRATION  

	%  $Revision$
 	%  was created 10-Nov-2021 14:12:37 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
 	%% It was developed on Matlab 9.11.0.1769968 (R2021b) for MACI64.  Copyright 2021 John Joowon Lee.
 	
	properties
 		
 	end

	methods 
		  
 		function this = Registration(varargin)
 			%% REGISTRATION
 			%  @param .

 			
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

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

