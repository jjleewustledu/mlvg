classdef Bids < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable  
	%% BIDS  

	%  $Revision$
 	%  was created 13-Nov-2021 14:58:16 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlpipeline/src/+mlpipeline.
 	%% It was developed on Matlab 9.11.0.1769968 (R2021b) for MACI64.  Copyright 2021 John Joowon Lee.
 	
	properties (Dependent)
        anatPath
        projPath
        derivativesPath
        destinationPath 		
        mriPath
        petPath
        sourcedataPath
        sourceAnatPath
        sourcePetPath
        subFolder
 	end

	methods 

        %% GET

        function g = get.anatPath(this)
            g = fullfile(this.derivativesPath, this.subFolder, 'anat', '');
        end
        function g = get.projPath(this)
            g = this.projPath_;
        end
        function g = get.derivativesPath(this)
            g = fullfile(this.projPath, 'derivatives', '');
        end
        function g = get.destinationPath(this)
            g = this.destPath_;
        end
        function g = get.mriPath(this)
            g = fullfile(this.derivativesPath, this.subFolder, 'mri', '');
        end
        function g = get.petPath(this)
            g = fullfile(this.derivativesPath, this.subFolder, 'pet', '');
        end
        function g = get.sourcedataPath(this)
            g = fullfile(this.projPath, 'sourcedata', '');
        end
        function g = get.sourceAnatPath(this)
            g = fullfile(this.sourcedataPath, this.subFolder, 'anat', '');
        end
        function g = get.sourcePetPath(this)
            g = fullfile(this.sourcedataPath, this.subFolder, 'pet', '');
        end
        function g = get.subFolder(this)
            g = this.subFolder_;
        end

        %%
		  
 		function this = Bids(varargin)
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

