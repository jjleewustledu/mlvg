classdef Reregistration < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
	%% REREGISTRATION  

	%  $Revision$
 	%  was created 07-Jun-2021 15:06:39 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
 	%% It was developed on Matlab 9.10.0.1669831 (R2021a) Update 2 for MACI64.  Copyright 2021 John Joowon Lee.
 	
	properties
 		radius = 1.1; % cf. sampleActivity()
 	end

	methods 
		  
 		function this = Reregistration(varargin)
 			%% REREGISTRATION  
            
            this.recentRewards_ = zeros(1, 3);
        end
        
        function [tform,centerlineOnTarget,reward] = pcregistermax(this, varargin)
            %  @param required tform is the initial tform, as rigid3d object.
            %  @param required centerline is a pointCloud, as numeric mlvg.Fung2013.U x 3.
            %  @param required target is a pointCloud.
            %  @return updated tform, as rigid3d object.
            %  @return centerlineOnTarget, as pointCloud.
            %  @return reward is the scalar optimized maximal activity sampled over the centerline.
            
            ip = inputParser;
            addRequired(ip, 'tform', @(x) isa(x, 'rigid3d'))
            addRequired(ip, 'centerline', @(x) isa(x, 'pointCloud'))
            addRequired(ip, 'target', @(x) isa(x, 'pointCloud'))
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            % init
            tform = ipr.tform;
            centerlineOnTarget = ipr.centerline;
            this.recentRewards_ = this.sampleActivity(ipr.centerline, ipr.target);
            idx = 1;
            while this.progressing()
                idx = idx + 1;
                [centerlineOnTarget,this.recentRewards_(idx)] = ...
                    this.pctransform();
            end
        end
        function [centerlineOnTarget,reward] = pctransform(this, tform, centerline, target)
            
            
        end
        function tf = progressing(this)
            
        end
        function a = sampleActivity(this, centerline, target)            
            %  @param required centerline is a pointCloud, as numeric mlvg.Fung2013.U x 3.
            %  @param required target is a pointCloud.
            %  @return a is activity from target sampled over the centerline.
            
            locs = centerline.Location;
            accum = [];
            for il = 1:size(locs, 1)
                idx = findNeighborsInRadius(target, centerline.Location(il,:), this.radius);
                selected = select(target, idx);
                accum(il) = mean(selected.Intensity); %#ok<AGROW>
            end
            a = mean(accum);
        end
    end 
    
    %% PRIVATE
    
    properties (Access = private)
        recentRewards_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end
