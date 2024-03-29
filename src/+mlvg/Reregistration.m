classdef Reregistration < handle
	%% REREGISTRATION 

	%  $Revision$
 	%  was created 07-Jun-2021 15:06:39 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
 	%% It was developed on Matlab 9.10.0.1669831 (R2021a) Update 2 for MACI64.  Copyright 2021 John Joowon Lee.
 	
	properties
        epsilon = 1e-2
        idxmax = 7 % max # of rigid body searches; cf. pcregistermax()
        mmppix
        Niterations = 3 % max # of translational searches
 		radius % cf. sampleActivity()
        target_bb % bounding box
        target_ic % imaging context
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

    properties (Dependent)
        ic_centerline1
        pc_centerline1
        stridea
        stridex
        tform1
        rewards1
    end

	methods 

        %% GET

        function g = get.ic_centerline1(this)
            ic_centerline1__ = copy(this.target_ic);
            ic_centerline1__.relocateToDerivativesFolder();
            ic_centerline1__.setPointCloud(this.centerline1_);            
            ic_centerline1__.fileprefix = strcat(this.target_ic.fileprefix, '_', clientname(true, 2));
            g = ic_centerline1__;
        end
        function g = get.pc_centerline1(this)
            g = this.centerline1_;
        end
        function g = get.stridea(~)
            g = 7.5; % 10 degrees may be too large
        end
        function g = get.stridex(this)
            g = 0.75/mean(this.mmppix); % 5 mm may be too large
        end
        function g = get.tform1(this)
            g = this.tform1_;
        end
        function g = get.rewards1(this)
            g = this.rewards1_;
        end

        %%

 		function this = Reregistration(varargin)
            ip = inputParser;
            addRequired(ip, 'target_ic');
            addParameter(ip, 'bounding_box', [], @(x) isa(x, 'mlaif.ArterialBoundingBox'))
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.target_ic = mlfourd.ImagingContext2(ipr.target_ic);
            this.target_ic.relocateToDerivativesFolder();
            this.target_bb = ipr.bounding_box;

            this.mmppix = this.target_ic.nifti.mmppix;
            this.radius = sqrt(2)*max(this.mmppix);
            this.voxelVolume = prod(this.mmppix)/1000;
        end
        
        function tuple = locationOfMaxIntensity(~, pc)
            %% returns [x y z] of pointCloud.Location that contains max pointCloud.Intensity.

            [~,loc] = max(pc.Intensity);
            tuple = pc.Location(loc,:);
        end
        function [tform1,centerline1,rewards1] = pcregistermax(this, varargin)
            %  @param required tform is the initial tform, as rigid3d|affine3d object.
            %  @param required centerline is a pointCloud, as numeric mlvg.Fung2013.U x 3.
            %  @param required target is a pointCloud.
            %  @return updated tform, as rigid3d object.
            %  @return centerline, as pointCloud.
            %  @return rewards is the array of maximal activities sampled over the centerline.
            
            ip = inputParser;
            addRequired(ip, 'tform', @(x) isa(x, 'rigid3d') || isa(x, 'affine3d'))
            addRequired(ip, 'centerline', @(x) isa(x, 'pointCloud'))
            addRequired(ip, 'target', @(x) isa(x, 'pointCloud'))
            parse(ip, varargin{:})
            ipr = ip.Results;
            tform = ipr.tform;
            this.centerline0_ = ipr.centerline;

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

            for iter = 1:this.Niterations
                % translational search
                idx_rew = idx_rew + 1;
                [tform,centerline_,rewards(idx_rew)] = this.pctransformmax3(tform, centerline_, target_, 1);
                disp(this)
                disp(rewards)
                fprintf('size(centerline.Location): %s\n', mat2str(size(centerline_.Location)))
            end
            for idx = 1:this.idxmax
                % rigid body search
                idx_rew = idx_rew + 1;
                [tform,centerline_,rewards(idx_rew)] = this.pctransformmax6(tform, centerline_, target_, idx);
                disp(this)
                disp(rewards)
                fprintf('size(centerline.Location): %s\n', mat2str(size(centerline_.Location)))
                if ~this.progressing(rewards)
                    break
                end
            end

            % reset center of native coords; update tform; target is not to be returned
            centerline1 = pctransform(centerline_, cform);
            tform1 = affine3d(cform.T * tform.T);
            rewards1 = rewards;
            
            this.centerline1_ = centerline1;
            this.tform1_ = tform1;
            this.rewards1_ = rewards1;
        end
        function [tform1,centerline1,rewards1] = pctransformmax6(this, tform, centerline, target, iteration)
            %% implements Fung & Carson, sec. 2.4, para. 3.:  5 samples per 6 affine d.o.f. ~ 15625 samples ~ 400 sec.

            tic
            centerline_ = copy(centerline); % don't clobber handle object
            da = 1/2^(iteration - 1);
            dx = da;
            stridea_ = this.stridea*da;
            Th = 2*stridea_;
            stridex_ = this.stridex*dx;
            L = 2*stridex_;
            for a1 = -Th:stridea_:Th % degrees
                for a2 = -Th:stridea_:Th
                    for a3 = -Th:stridea_:Th
                        for x1 = -L:stridex_:L % voxel widths
                            for x2 = -L:stridex_:L
                                for x3 = -L:stridex_:L

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
                                    X1 = eye(4);
                                    X1(4,1:3) = [x1 x2 x3];
                                    tform_ = affine3d(X1 * A3 * A2 * A1);
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

            tform1 = affine3d(this.tform_star.T * tform.T);
            centerline1 = this.centerline_star;
            rewards1 = this.sa_star;
            fprintf('mlvg.Reregistration.pctransformmax6: ')
            toc
        end
        function [tform1,centerline1,rewards1] = pctransformmax3(this, tform, centerline, target, iteration)
            %% implements Fung & Carson, sec. 2.4, para. 3.:  9 samples per 3 affine d.o.f. ~ 729 samples ~ 15 sec.

            tic
            centerline_ = copy(centerline); % don't clobber handle object
            dx = 1/2^(iteration - 1);
            stridex_ = this.stridex*dx;
            L = 4*stridex_;
            for x1 = -L:stridex_:L % voxel widths
                for x2 = -L:stridex_:L
                    for x3 = -L:stridex_:L

                        X1 = eye(4);
                        X1(4,1:3) = [x1 x2 x3];
                        tform_ = affine3d(X1);
                        centerline__ = pctransform(centerline_, tform_); % centerline trial
                        sa_ = this.sampleActivity(centerline__, target);
                        if sa_ > this.sa_star % update starred objects
                            this.tform_star = tform_;
                            this.centerline_star = copy(centerline__);
                            this.sa_star = sa_;
                            this.a1_star = 0;
                            this.a2_star = 0;
                            this.a3_star = 0;
                            this.x1_star = x1;
                            this.x2_star = x2;
                            this.x3_star = x3;
                        end
                    end
                end
            end

            tform1 = affine3d(this.tform_star.T * tform.T);
            centerline1 = this.centerline_star;
            rewards1 = this.sa_star;
            fprintf('mlvg.Reregistration.pctransformmax3: ')
            toc
        end        
        function h = pcshow(this, varargin)
            h = figure;
            hold on; 
            pcshow(this.target_ic.pointCloud());
            pcshow(this.centerline1_.Location, 'm', 'MarkerSize', 3);
            pcshow(this.centerline0_.Location, '+g', 'MarkerSize', 3);
            hold off;
            view(-210, 15) % tuned for orientstd
            daspect([1 1 1]);

            %set(h, 'InvertHardCopy', 'off');
            %set(h, 'Color', [0 0 0]);
            %set(gca, 'Color', [0 0 0]);
            xlabel('x (mm)', 'FontSize', 14)
            ylabel('y (mm)', 'FontSize', 14)
            zlabel('z (mm)', 'FontSize', 14)
            set(gca, 'xcolor', [1 1 1])
            set(gca, 'ycolor', [1 1 1])
            set(gca, 'zcolor', [1 1 1])

            title(strcat(clientname(true, 2)), 'Interpreter', 'none', 'FontSize', 14);
        end
        function tf = progressing(this, rewards)
            %% progression is reward increasing by at least an epsilon of the previous reward

            if length(rewards) < this.Niterations + 1
                tf = true;
                return
            end
            tf = rewards(end) - rewards(end-1) > this.epsilon*rewards(end-1);
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

    %% PROTECTED

    properties (Access = protected)
        centerline0_
        centerline1_
        tform1_
        rewards1_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end
