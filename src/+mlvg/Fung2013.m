classdef Fung2013 < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable  
    %% FUNG2013 implements
    %  Edward K Fung and Richard E Carson.  Cerebral blood flow with [15O]water PET studies using 
    %  an image-derived input function and MR-defined carotid centerlines.  
    %  Phys. Med. Biol. 58 (2013) 1903–1923.  doi:10.1088/0031-9155/58/6/1903
    
    %  $Revision$
 	%  was created 22-Mar-2021 22:11:00 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
 	%% It was developed on Matlab 9.9.0.1592791 (R2020b) Update 5 for MACI64.  Copyright 2021 John Joowon Lee.
    
    properties
        BBBuf % extra voxels padded to coords to create convex bounding box
        T1w_ic
        centerlines_ics % L, R in cell
        centerlines_pcs % L, R in cell
        coords % 4 points at corners
        coords_bb % {x1:xN, y1:yN, z1:zN} for bounding box
        corners_ic
        cornersb_ic
        plotmore % show more plots for QA
        ploton % show final results
        registration % struct
            % tform
            % centerlineOnTarget
            % rmse 
            % targets are averages of frames containing 10-25 pcnt, 10-50 pcnt, 10-75 pcnt of max emissions
        reregistration % struct
            % tform
            % centerlineOnTarget
            % reward 
            % targets are averages of frames containing 10-25 pcnt, 10-50 pcnt, 10-75 pcnt of max emissions
        segmentation_ic % contains solid 3D volumes for carotids
        wmparc_ic
        
        % for B-splines in mlvg.Hunyadi2021
        k = 4
        t = [0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
        %t = [0 0 0 0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1 1 1]
        U % # samples for bspline_deboor()
        Cs % curves in cell, 3 x this.U
        Ps % control points in cell, 3 x M
    end
    
    properties (Dependent)
        anatPath
        CCIRPath
        derivativesPath
        Nx
        Ny
        Nz
        mriPath
        petPath
        sourcedataPath
        subFolder
    end

    methods
        
        %% GET
        
        function g = get.anatPath(this)
            g = fullfile(this.sourcedataPath, this.subFolder, 'anat', '');
        end
        function g = get.CCIRPath(this)
            g = fileparts(fileparts(fileparts(this.petPath_)));
        end
        function g = get.derivativesPath(this)
            g = fullfile(this.CCIRPath, 'derivatives', '');
        end
        function g = get.Nx(this)
            g = size(this.T1w_ic, 1);
        end
        function g = get.Ny(this)
            g = size(this.T1w_ic, 2);
        end
        function g = get.Nz(this)
            g = size(this.T1w_ic, 3);
        end
        function g = get.mriPath(this)
            g = fullfile(this.derivativesPath, this.subFolder, 'mri', '');
        end
        function g = get.petPath(this)
            g = this.petPath_;
        end
        function g = get.sourcedataPath(this)
            g = fullfile(this.CCIRPath, 'sourcedata', '');
        end
        function g = get.subFolder(this)
            g = fileparts(this.petPath_);
            g = basename(g);
        end
        
        %%
        
        function this = Fung2013(varargin)
            %% FUNG2013
            %  @param coords from fsleyes [ x y z; ... ], [ [RS]; [LS]; [RI]; [LI] ].
            %  @param iterations ~ 130.
            %  @param smoothFactor ~ 0.
            
            ip = inputParser;
            addOptional(ip, 'petPath', pwd, @isfolder)
            addParameter(ip, 'ploton', true, @islogical)
            addParameter(ip, 'plotmore', false, @islogical)
            addParameter(ip, 't1w', 'sub-*_T1w.nii.gz', @ischar)
            addParameter(ip, 'coords', [], @isnumeric)
            addParameter(ip, 'BBBuf', [10 10 4], @isnumeric)
            addParameter(ip, 'iterations', 100, @isscalar)
            addParameter(ip, 'smoothFactor', 0, @isscalar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.petPath_ = ipr.petPath;
            this.ploton = ipr.ploton;
            this.plotmore = ipr.plotmore;
            this.coords = ipr.coords;
            this.BBBuf = ipr.BBBuf;
            
            % gather requirements
            t1w = globT(fullfile(this.anatPath, ipr.t1w));
            assert(isfile(t1w{1}))
            this.T1w_ic = mlfourd.ImagingContext2(t1w{1});
            this.hunyadi_ = mlvg.Hunyadi2021();
            
            % build segmentation
            this.buildCorners(this.coords);
            this.buildSegmentation(ipr.iterations, 'smoothFactor', ipr.smoothFactor);
        end
        function this = buildCorners(this, varargin)
            %  @param coords is [x y z; x2 y2 z2; x3 y3 z3; x4 y4 z4] | empty.
            %         coords is [ [RS]; [LS]; [RI]; [LI] ].
            %  @return this.{corners*_ic, *_range}
            
            %  f = mlvg.Fung2013
            %  f.buildCorners([158 122 85; 96 126 88; 156 116 27; 101 113 28])
            %  158, 122, 85
            %  96, 126, 88
            %  156, 116, 27
            %  101, 113, 28

            ip = inputParser;
            addOptional(ip, 'coords', this.coords, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.coords = ipr.coords;
            cc = num2cell(this.coords);
            
            % pick corners
            if isempty(this.coords)  
                disp('No coords for carotids are available.  Please find coords in the T1w and provide to the constructor.')
                assert(~isempty(this.T1w_ic), 'Oops:  No T1w is available.  Please provide information for T1w to the constructor.')
                this.T1w_ic.fsleyes
                error('mlvg:Fung2013', ...
                    'No coords for carotids were available.  Please provide carotid coords to the constructor.')
            else                
                assert(all(size(this.coords) == [4 3]))
            end
            
            % build ImagingContexts with corners, and also with blurring, binarizing
            this.corners_ic = this.T1w_ic.zeros;
            this.corners_ic.fileprefix = 'corners_on_T1w';
            nii = this.corners_ic.nifti;
            nii.img(cc{1,:}) = 1;
            nii.img(cc{2,:}) = 1;
            nii.img(cc{3,:}) = 1;
            nii.img(cc{4,:}) = 1;
            this.corners_ic = mlfourd.ImagingContext2(nii);
            assert(4 == dipsum(this.corners_ic))
            
            this.cornersb_ic = this.corners_ic.blurred(1);
            this.cornersb_ic = this.cornersb_ic.numgt(0.001);
            this.cornersb_ic.fileprefix = 'corners_on_T1w_spheres';
            
            % build coords_bb
            sz = size(this.T1w_ic);
            for m = 1:3
                this.coords_bb{m} = (min(ipr.coords(:,m)) - this.BBBuf(m)):(max(this.coords(:,m)) + this.BBBuf(m) + 1);
                bb = this.coords_bb{m};
                this.coords_bb{m} = bb(bb >=1 & bb <= sz(m));
            end
        end
        function this = buildSegmentation(this, varargin)
            %  @param optional iterations ~ 100.
            %  @param smoothFactor ~ 0.
            %  @return this.segmentation_ic.
            
            ip = inputParser;
            addOptional(ip, 'iterations', 100, @isscalar)
            addParameter(ip, 'smoothFactor', 0, @isscalar)
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            T1wb_img = this.T1w_ic.nifti.img(this.coords_bb{:});
            cornersb_img = this.cornersb_ic.nifti.img(this.coords_bb{:});
            
            % call snakes, viz., iterate
            ac = activecontour(T1wb_img, cornersb_img, ipr.iterations, 'Chan-Vese', 'SmoothFactor', ipr.smoothFactor);
            if this.plotmore
                this.plotSegmentation(ac)
                title(sprintf('iterations %i, smooth %g', ipr.iterations, ipr.smoothFactor))
            end

            % fit back into T1w
            ic = this.T1w_ic.zeros; 
            ic.fileprefix = 'segmentation_on_T1w';
            nii = ic.nifti;
            nii.img(this.coords_bb{:}) = ac;
            nii.save()
            this.segmentation_ic = mlfourd.ImagingContext2(nii);
        end
        function this = buildCenterlines(this)
            img = logical(this.segmentation_ic);
            imgL = img(1:ceil(this.Nx/2),:,:);
            imgR = zeros(size(img));
            imgR(ceil(this.Nx/2)+1:end,:,:) = img(ceil(this.Nx/2)+1:end,:,:);
            [pcL,CL,PL] = this.buildCenterline(imgL, '_left');
            [pcR,CR,PR] = this.buildCenterline(imgR, '_right');
            this.centerlines_pcs = {pcL pcR};
            this.Cs = {CL CR};
            this.Ps = {PL PR};
        end
        function [pc,C,P] = buildCenterline(this, img, tag)
            assert(ischar(tag))
            idx = find(img);
            [X,Y,Z] = ind2sub(size(img), idx);             
            M(1,:) = X'; % M are ints cast as double
            M(2,:) = Y';
            M(3,:) = Z';
            this.U = 2*length(Z);
            
            P = bspline_estimate(this.k, this.t, M); % double
            C = bspline_deboor(this.k, this.t, P, this.U); % double, ~2x oversampling for Z
            pc = pointCloud(C');
            
            if this.plotmore
                figure;
                hold all;
                plot3(M(1,:), M(2,:), M(3,:), 'k.');
                plot3(P(1,:), P(2,:), P(3,:), 'b');
                plot3(C(1,:), C(2,:), C(3,:), 'r');
                legend('segmentation', 'control points', 'curve', ...
                    'Location', 'Best');
                hold off;
                figure;
                pcshow(pc)
            end
        end
        function this = buildRegistrationTargets(this, dyn_ic)
            this.wmparc_ic = mlfourd.ImagingContext2(fullfile(this.mriPath, 'wmparc_on_T1w.nii.gz'));
            dyn_avgxyz = dyn_ic.volumeAveraged(logical(this.wmparc_ic));
            dyn_max = dipmax(dyn_avgxyz);
            img = dyn_avgxyz.nifti.img;
            [~,it10] = max(img > 0.1*dyn_max);
            [~,it25] = max(img > 0.25*dyn_max);
            [~,it50] = max(img > 0.5*dyn_max);
            [~,it75] = max(img > 0.75*dyn_max);
            
            this.registration.target_ics{1} = dyn_ic.timeAveraged(it10:it25);
            this.registration.target_ics{2} = dyn_ic.timeAveraged(it10:it50);
            this.registration.target_ics{3} = dyn_ic.timeAveraged(it10:it75);
            for i = 1:3
                this.registration.target_ics{i} = this.registration.target_ics{i}.blurred(2.5);
            end
        end
        function [ics,dyn_ics] = call(this, varargin)
            %% CALL
            %  @param optional toglob, e.g., 'sub-*Dynamic*_on_T1w.nii.gz'
            %  @return {ImagingContext2 objects for centerlines}
            %  @return {ImagingContext2 objects for globbed}
            
            ip = inputParser;
            addOptional(ip, 'toglob', 'sub-*Dynamic*_on_T1w.nii.gz', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            % build centerlines
            this.buildCenterlines()

            ics = {};
            dyn_ics = {};
            for niis = globT(ipr.toglob)                
                dyn_ic = mlfourd.ImagingContext2(niis{1});                
                this.buildRegistrationTargets(dyn_ic)
                this.registerCenterlines('alg', 'cpd')
                ic = this.pointCloudsToIC();
                ic.filepath = dyn_ic.filepath;
                ic.fileprefix = [dyn_ic.fileprefix '_idifmask'];
                ic.save()
                ics = [ics {ic}]; %#ok<AGROW>
                dyn_ics = [dyn_ics {dyn_ic}]; %#ok<AGROW>
            end
        end
        function h = plotRegistered(this, varargin)
            % @param required target, pointCloud.
            % @param required centerlineOnTarget, pointCloud.
            % @param required centerline, pointCloud.
            % @param required larerality, in {'' 'l' 'L' 'r' 'R'}.
            
            ip = inputParser;
            addRequired(ip, 'target', @(x) isa(x, 'pointCloud'))
            addRequired(ip, 'centerlineOnTarget', @(x) isa(x, 'pointCloud'))
            addRequired(ip, 'centerline', @(x) isa(x, 'pointCloud'))
            addRequired(ip, 'laterality', @(x) ismember(x, {'', 'l', 'L', 'r', 'R'}))
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            if this.plotmore
                %figure;
                %pcshow(ipr.centerline)
                %title(sprintf('centerline %s', upper(ipr.laterality)))
                figure
                pcshow(ipr.centerlineOnTarget)
                title(sprintf('centerlineOnTarget %s', upper(ipr.laterality)))
            end
            h = figure;
            pcshowpair(ipr.target, ipr.centerlineOnTarget, 'VerticalAxis', 'Z')
            title(sprintf('target (magenta) + centerlineOnTarget (green) %s', upper(ipr.laterality)))
        end
        function h = plotSegmentation(~, ac)
            %  @param required activecontour result.
            
            h = figure;
            p = patch(isosurface(double(ac)));
            p.FaceColor = 'red';
            p.EdgeColor = 'none';
            daspect([1 1 1])
            camlight;
            lighting phong            
        end
        function ic = pointCloudsToIC(this, varargin)
                icL = this.pointCloudToIC(this.registration.centerlineOnTarget{1}, varargin{:});
                icL = icL.imdilate(strel('sphere', 2));
                icR = this.pointCloudToIC(this.registration.centerlineOnTarget{2}, varargin{:});
                icR = icR.imdilate(strel('sphere', 2));
                ic = icL + icR;
                assert(1 == dipmax(ic))
        end
        function ic = pointCloudToIC(this, pc, varargin)
            ip = inputParser;
            addRequired(ip, 'pc', @(x) isa(x, 'pointCloud'))
            addOptional(ip, 'fileprefix', 'pointClouudToIC', @ischar)
            parse(ip, pc, varargin{:})
            ipr = ip.Results;
            
            ic = this.T1w_ic.zeros();
            ifc = ic.nifti;
            ifc.fileprefix = ipr.fileprefix;
            X = round(pc.Location(:,1));
            Y = round(pc.Location(:,2));
            Z = round(pc.Location(:,3));
            ind = sub2ind(size(ifc), X, Y, Z);
            ifc.img(ind) = 1;
            ic = mlfourd.ImagingContext2(ifc);
        end
        function this = registerCenterlines(this, varargin)
            %  @param thresh applies to ic3d.  Default is 25000.
            %  @param alg is from {'ndt', 'icp', 'cpd'}.
            %  @param gridStep preprocesses pcregister* methods.
            
            assert(~isempty(this.centerlines_pcs))
            this.registerCenterline(  this.centerlines_pcs{1}, varargin{:}, 'laterality', 'L') 
            this.reregisterCenterline(this.centerlines_pcs{1}, varargin{:}, 'laterality', 'L')            
            this.registerCenterline(  this.centerlines_pcs{2}, varargin{:}, 'laterality', 'R')           
            this.reregisterCenterline(this.centerlines_pcs{2}, varargin{:}, 'laterality', 'R')
        end
        function this = registerCenterline(this, centerline, varargin)
            %  @param required centerline is a pointCloud.
            %  @param optional ic3d is an ImagingContext2 for PET averaged over early times for bolus arrival.
            %         Default is this.registration.target_ics{2}.
            %  @param thresh applies to ic3d.  Default is 25000.
            %  @param alg is from {'ndt', 'icp', 'cpd'}.
            %  @param gridStep preprocesses pcregister* methods.
            %  @param laterality is in {'R' 'L'}.
            
            ip = inputParser;
            addRequired(ip, 'centerline', @(x) isa(x, 'pointCloud'))
            addOptional(ip, 'ic3d', this.registration.target_ics{2}, @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'thresh', 25000, @isnumeric)
            addParameter(ip, 'alg', 'cpd', @(x) ismember(x, {'ndt', 'icp', 'cpd'}))
            addParameter(ip, 'gridStep', 1, @isscalar)
            addParameter(ip, 'laterality', '', @(x) ismember(x, {'', 'l', 'L', 'r', 'R'})) % L has indices < Nx/2
            parse(ip, centerline, varargin{:})
            ipr = ip.Results;
            ipr.ic3d = this.maskInBoundingBox(ipr.ic3d, ipr.laterality);
            if isempty(ipr.thresh)
                ipr.thresh = 0.1*dipmax(ipr.ic3d);
            end
            target = pointCloud(ipr.ic3d, 'thresh', ipr.thresh); 
            if this.ploton
                figure
                pcshow(target)
                title(sprintf('target %s', upper(ipr.laterality)))           
            end
            
            switch ipr.alg
                case 'ndt'
                    [tform,centerlineOnTarget,rmse] = pcregisterndt(ipr.centerline, target, ipr.gridStep, ...
                        'Tolerance', [0.01 0.05]);
                case 'icp'
                    if ipr.gridStep ~= 1
                        ipr.centerline = pcdownsample(ipr.centerline, 'gridAverage', ipr.gridStep);
                    end
                    [tform,centerlineOnTarget,rmse] = pcregistericp(ipr.centerline, target, ...
                        'Extrapolate', true, 'Tolerance', [0.01 0.01]);
                case 'cpd'
                    if ipr.gridStep ~= 1
                        ipr.centerline = pcdownsample(ipr.centerline, 'gridAverage', ipr.gridStep);
                    end
                    [tform,centerlineOnTarget,rmse] = pcregistercpd(ipr.centerline, target, ...
                        'Transform', 'Rigid', 'MaxIterations', 100, 'Tolerance', 1e-7); % 'InteractionSigma', 2
                otherwise
                    error('mlvg:ValueError', ...
                        'Fung2013.registerCenterlines.ipr.alg == %s', ipr.alg)
            end
            
            idx = strcmpi(ipr.laterality, 'R') + 1; % idx == 1 <-> left
            this.registration.tform{idx} = tform;
            this.registration.centerlineOnTarget{idx} = centerlineOnTarget;
            this.registration.rmse{idx} = rmse;
       
            if this.ploton
                this.plotRegistered(target, centerlineOnTarget, ipr.centerline, ipr.laterality)
            end
        end
        function this = reregisterCenterline(this, varargin)
            %  @param required centerline is pointCloud, this.U x 3.
            %  @param optional ic3d is an ImagingContext2 for PET averaged over early times for bolus arrival.
            %         Default is this.registration.target_ics{2}.
            %  @param thresh applies to ic3d.  Default is 25000.
            %  @param laterality is in {'R' 'L'}.
            
            ip = inputParser;
            addRequired(ip, 'centerline', @(x) isa(x, 'pointCloud'))          
            addOptional(ip, 'ic3d', this.registration.target_ics{2}, @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'thresh', 25000, @isnumeric)
            addParameter(ip, 'laterality', '', @(x) ismember(x, {'', 'l', 'L', 'r', 'R'})) % L has indices < Nx/2
            parse(ip, centerline, varargin{:})
            ipr = ip.Results;
            ipr.ic3d = this.maskInBoundingBox(ipr.ic3d, ipr.laterality);
            if isempty(ipr.thresh)
                ipr.thresh = 0.1*dipmax(ipr.ic3d);
            end
            target = pointCloud(ipr.ic3d, 'thresh', ipr.thresh);
            if this.ploton
                figure
                pcshow(target)
                title(sprintf('target %s', upper(ipr.laterality)))           
            end            
            
            idx = strcmpi(ipr.laterality, 'R') + 1; % idx == 1 <-> left
            rr = mlvg.Reregistration();
            [tform,centerlineOnTarget,reward] = rr.pcregistermax( ...
                this.registration.tform{idx}, ipr.centerline, target);
            this.reregistration.tform{idx} = tform;
            this.reregistration.centerlineOnTarget{idx} = centerlineOnTarget;
            this.reregistration.reward{idx} = reward;
            
            if this.ploton
                this.plotRegistered( ...
                    target, centerlineOnTarget, ipr.centerline, ipr.laterality)
            end
        end
    end
    
    %% PROTECTED    
    
    methods (Access = protected)
        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
            that.hunyadi_ = copy(this.hunyadi_);
        end
    end
    
    %% PRIVATE
    
    properties (Access = private)
        hunyadi_
        petPath_
    end
    
    methods (Access = private)
        function ic = maskInBoundingBox(this, ic, laterality)
            ifc = ic.nifti;
            img = zeros(size(ifc));
            bb = this.coords_bb;
            img(bb{1},bb{2},bb{3}) = ifc.img(bb{1},bb{2},bb{3});
            if ~isempty(laterality)
                if strcmpi(laterality, 'L') % L has indices < Nx/2
                    img(ceil(this.Nx/2)+1:end,:,:) = zeros(ceil(this.Nx/2),this.Ny,this.Nz);
                else
                    img(1:ceil(this.Nx/2),:,:) = zeros(ceil(this.Nx/2),this.Ny,this.Nz);                
                end
            end
            ifc.img = img;
            ifc.fileprefix = [ifc.fileprefix '_maskInBoundingBox'];
            ic = mlfourd.ImagingContext2(ifc);
        end
    end
end