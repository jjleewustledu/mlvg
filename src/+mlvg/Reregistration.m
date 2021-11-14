classdef Reregistration < handle
	%% REREGISTRATION 

	%  $Revision$
 	%  was created 07-Jun-2021 15:06:39 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
 	%% It was developed on Matlab 9.10.0.1669831 (R2021a) Update 2 for MACI64.  Copyright 2021 John Joowon Lee.
 	
	properties
        epsilon = 1e-2
        idxmax = 8 % cf. pcregistermax()
        mmppix
        Niterations = 1 % # of iterated searches per specified search grid (dangle, dpos)
 		radius % cf. sampleActivity()
        target_ic
        voxelVolume % mL

        a1_star = [];
        a2_star = [];
        a3_star = [];
        x1_star = [];
        x2_star = [];
        x3_star = [];            
        tform_star = [];
        centerline_star = [];
        sa_star = 0;
 	end

	methods 
		  
 		function this = Reregistration(varargin)
            ip = inputParser;
            addRequired(ip, 'target_ic', @(x) isa(x, 'mlfourd.ImagingContext2'))
            parse(ip, varargin{:})
            ipr = ip.Results;

            this.target_ic = ipr.target_ic;
            this.mmppix = this.target_ic.nifti.mmppix;
            this.radius = min(this.mmppix) + this.epsilon;
            this.voxelVolume = prod(this.mmppix)/1000;
        end
        
        function tuple = locationOfMaxIntensity(~, pc)
            %% returns [x y z] of pointCloud.Location that contains max pointCloud.Intensity.

            [~,loc] = max(pc.Intensity);
            tuple = pc.Location(loc,:);
        end
        function [tform,centerline1,rewards] = pcregistermax(this, varargin)
            %  @param required tform is the initial tform, as rigid3d object.
            %  @param required centerline is a pointCloud, as numeric mlvg.Fung2013.U x 3.
            %  @param required target is a pointCloud.
            %  @return updated tform, as rigid3d object.
            %  @return centerline, as pointCloud.
            %  @return rewards is the array of maximal activities sampled over the centerline.
            
            ip = inputParser;
            addRequired(ip, 'tform', @(x) isa(x, 'rigid3d'))
            addRequired(ip, 'centerline', @(x) isa(x, 'pointCloud'))
            addRequired(ip, 'target', @(x) isa(x, 'pointCloud'))
            parse(ip, varargin{:})
            ipr = ip.Results;
            tform = ipr.tform;

            % set center of transformations to be the coord of max sampleActivity()
            C = eye(4);
            C(4,1:3) = this.locationOfMaxIntensity(ipr.target); % coord tuple of hotspot
            cform = affine3d(C);
            cform_inv = cform.invert;
            tform = affine3d(cform_inv.T * tform.T);
            
            % implement pctransformmax using sampleActivity() as objective
            centerline_ = pctransform(ipr.centerline,  cform_inv);
            target_ = pctransform(ipr.target,  cform_inv);
            rewards = this.sampleActivity(centerline_, target_);
            idx_rew = 1;
            for idx = 1:this.idxmax
                for iter = 1:this.Niterations
                    idx_rew = idx_rew + 1;
                    [tform,centerline_,rewards(idx_rew)] = this.pctransformmax(tform, centerline_, target_, idx);
                    disp(this)
                    disp(rewards)
                end
                if ~this.progressing(rewards)
                    break
                end
            end

            % reset center of native coords; update tform; target is not to be returned
            centerline1 = pctransform(centerline_, cform);
            tform = affine3d(cform.T * tform.T);
        end
        function [tform,centerline1,reward] = pctransformmax(this, tform, centerline, target, iteration)
            %% implements Fung & Carson, sec. 2.4, para. 3., but only 5 samples per 6 affine d.o.f. ~ 15625 samples.

            tic
            centerline_ = copy(centerline); % don't clobber handle object
            da = 1/2^(iteration - 1);
            dx = da;
            for a1 = -20*da:10*da:20*da % degrees
                for a2 = -20*da:10*da:20*da
                    for a3 = -20*da:10*da:20*da
                        for x1 = -10*dx:5*dx:10*dx % voxel widths
                            for x2 = -10*dx:5*dx:10*dx
                                for x3 = -10*dx:5*dx:10*dx

                                    A1 = [1 0 0 0; ...
                                          0  cosd(a1) sind(a1) 0; ...
                                          0 -sind(a1) cosd(a1) 0; ...
                                          0 0 0 1];
                                    A2 = [cosd(a2) 0 -sind(a2) 0; ...
                                          0 1 0 0; ...
                                          sind(a2) 0  cosd(a2) 0;
                                          0 0 0 1];
                                    A3 = [ cosd(a3) sind(a3) 0 0; ...
                                          -sind(a3) cosd(a3) 0 0; ...
                                           0 0 1 0;
                                           0 0 0 1];
                                    X = eye(4);
                                    X(4,1:3) = [x1 x2 x3];
                                    tform_ = affine3d(X * A3 * A2 * A1);
                                    centerline__ = pctransform(centerline_, tform_); % centerline trial
                                    sa_ = this.sampleActivity(centerline__, target);
                                    if sa_ > this.sa_star % update starred objects
                                        this.tform_star = tform_;
                                        this.centerline_star = copy(centerline__);
                                        this.sa_star = sa_;
                                        this.a1_star = a1;
                                        this.a2_star = a2;
                                        this.a3_star = a3;
                                        this.x1_star = x1;
                                        this.x2_star = x2;
                                        this.x3_star = x3;
                                    end
                                end
                            end
                        end
                    end
                end
            end

            tform = affine3d(this.tform_star.T * tform.T);
            centerline1 = this.centerline_star;
            reward = this.sa_star;
            fprintf('mlvg.Reregistration.pctransformmax: ')
            toc
        end
        function tf = progressing(this, rewards)
            %% progression is reward increasing by at least an epsilon of the previous reward

            if length(rewards) < this.Niterations + 1
                tf = true;
                return
            end
            tf = all(rewards(end) - rewards(end-this.Niterations:end-1) > this.epsilon*rewards(end));
        end
        function a = sampleActivity(this, centerline, target)
            %  @param required centerline is a pointCloud, numeric of size mlvg.Fung2013.U x 3.
            %  @param required target is a pointCloud of time-averaged PET.
            %  @return a is the scalar activity from the target sampled by the centerline (Bq).
            
            locs = centerline.Location; % M x 3
            accum = zeros(1, size(locs, 1));
            R = this.radius;
            for il = 1:size(locs, 1)
                neighIndices = findNeighborsInRadius(target, locs(il,:), R);
                if ~isempty(neighIndices)
                    selected = select(target, neighIndices);
                    accum(il) = mean(selected.Intensity); % moving average of neighbors
                end
            end
            a = sum(accum) * this.voxelVolume;
        end
    end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end
