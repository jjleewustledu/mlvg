classdef VisionBids < handle & mlpipeline.Bids
	%% VISIONBIDS  

	%  $Revision$
 	%  was created 13-Nov-2021 14:57:47 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
 	%% It was developed on Matlab 9.11.0.1769968 (R2021b) for MACI64.  Copyright 2021 John Joowon Lee.

    properties
        toglob
    end

	methods 		  
 		function this = VisionBids(varargin)
            this = this@mlpipeline.Bids(varargin{:});
            this.toglob = fullfile(this.petPath, 'sub-*Dynamic*_on_T1w.nii.gz');
        end
    end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

