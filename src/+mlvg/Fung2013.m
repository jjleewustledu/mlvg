classdef Fung2013 < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable  
    %% FUNG2013 implements
    %  Edward K Fung and Richard E Carson.  Cerebral blood flow with [15O]water PET studies using 
    %  an image-derived input function and MR-defined carotid centerlines.  
    %  Phys. Med. Biol. 58 (2013) 1903–1923.  doi:10.1088/0031-9155/58/6/1903
    %  See also:  mlvg.Registration, mlvg.Reregistration
    
    %  $Revision$
 	%  was created 22-Mar-2021 22:11:00 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
 	%% It was developed on Matlab 9.9.0.1592791 (R2020b) Update 5 for MACI64.  Copyright 2021 John Joowon Lee.
    
    properties
        alg = 'fung' % prefer 'fung'; try 'cpd'
        BBBuf % extra voxels padded to coords to create convex bounding box
        centerlines_ics % L, R in cell
        centerlines_pcs % L, R in cell
        coords % 4 points at corners
        coords_bb % {x1:xN, y1:yN, z1:zN} for bounding box
        corners_ic
        cornersb_ic
        dilationRadius = 2 % Fung reported best results with radius ~ 2.5, but integers may be faster
        dyn_fileprefix % fileprefix for PET
        dyn_label % label for PET, useful for figures
        it10
        it25
        it50
        it75
        plotdebug % show debugging plots
        plotmore % show more plots for QA
        ploton % show final results
        registration % struct
            % tform
            % centerlineOnTarget
            % rmse 
            % target_ics are averages of frames containing 10-25 pcnt, 10-50 pcnt, 10-75 pcnt of max emissions
        segmentation_ic % contains solid 3D volumes for carotids
        T1w_ic
        taus % containers.Map
        times % containers.Map
        timesMid % containers.Map
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
        projPath
        derivativesPath
        destinationPath
        NCenterlineSamples % 1 voxel/mm for coarse representation of b-splines
        Nx
        Ny
        Nz
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
        function g = get.NCenterlineSamples(this)
            rngz = max(this.coords_bb{3}) - min(this.coords_bb{3});
            g = ceil(rngz/this.T1w_ic.nifti.mmppix(3)); % sample to 1 mm, or 1 voxels/mm
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
        
        function this = Fung2013(varargin)
            %% FUNG2013
            %  @param optional destPath is the path for writing outputs.  Default is pwd.  
            %         Must specify project ID & subject ID
            %  @param t1w is a filename string to glob.
            %  @param coords from fsleyes [ x y z; ... ], [ [RS]; [LS]; [RI]; [LI] ].
            %  @param iterations ~ 130.
            %  @param smoothFactor ~ 0.
            
            ip = inputParser;
            addOptional(ip, 'destPath', pwd, @isfolder)
            addParameter(ip, 'ploton', true, @islogical)
            addParameter(ip, 'plotmore', true, @islogical)
            addParameter(ip, 'plotdebug', false, @islogical)
            addParameter(ip, 't1w', 'sub-*_T1w.nii.gz', @ischar)
            addParameter(ip, 'coords', [], @isnumeric)
            addParameter(ip, 'BBBuf', [16 16 4], @isnumeric)
            addParameter(ip, 'iterations', 100, @isscalar)
            addParameter(ip, 'smoothFactor', 0, @isscalar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.parseDestinationPath(ipr.destPath);
            this.ploton = ipr.ploton;
            this.plotmore = ipr.plotmore;
            this.plotdebug = ipr.plotdebug;
            this.coords = ipr.coords;
            this.BBBuf = ipr.BBBuf;
            
            % gather requirements
            t1w = globT(fullfile(this.sourceAnatPath, ipr.t1w));
            assert(isfile(t1w{1}))
            this.T1w_ic = mlfourd.ImagingContext2(t1w{1});
            this.hunyadi_ = mlvg.Hunyadi2021();
            this.buildCorners(this.coords);

            % some timing objects
            this.taus = containers.Map;
            this.taus('CO') = [15 60 60 60 60 60];
            this.taus('HO') = [3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 5 5 5 5 5 5 10 10 10 10 10 10 10 10 30 30 30 30 30 30];
            this.taus('OO') = [3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 5 5 5 5 5 5 10 10 10 10 10 10 10 10 30 30 30 30 30 30];
            this.taus('FDG') = [5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 20 20 20 20 20 20 20 20 20 60 60 60 60 60 60 60 60 60 60 300 300 300 300 300 300 300 300 300];

            this.times = containers.Map;
            for key = this.taus.keys
                 this.times(key{1}) = [0 cumsum(this.taus(key{1}))];
            end

            this.timesMid = containers.Map;
            for key = this.taus.keys
                taus_ = this.taus(key{1});
                times_ = this.times(key{1});
                this.timesMid(key{1}) = times_(1:length(taus_)) + taus_/2;
            end
        end
        function this = buildCorners(this, varargin)
            %% BUILDCORNERS builds representations of the bounding box as images and coord ranges.
            %  As needed, it launches fsleyes for manual selection of bounding box corners.
            %  @param coords is [x y z; x2 y2 z2; x3 y3 z3; x4 y4 z4] | empty.
            %         coords is [ [RS]; [LS]; [RI]; [LI] ].
            %  @return this.corners*_ic, which represent corners of the bounding box with unit voxels in arrays of zeros.
            %  @return this.coords_bb, which are row arrays for bases [x y z] that describe the range of bounding box voxels.
            
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
            %% segments the arterial path using activecontour() with the 'Chan-Vese' method.            
            %  @param optional iterations ~ 100.
            %  @param smoothFactor ~ 0.
            %  @return this.segmentation_ic.
            
            ip = inputParser;
            addOptional(ip, 'iterations', 100, @isscalar)
            addParameter(ip, 'smoothFactor', 0, @isscalar)
            parse(ip, varargin{:})
            ipr = ip.Results;

            if ~isempty(this.segmentation_ic)
                return
            end
                        
            T1wb_img = this.T1w_ic.nifti.img(this.coords_bb{:});
            cornersb_img = this.cornersb_ic.nifti.img(this.coords_bb{:});
            
            % call snakes, viz., iterate
            ac = activecontour(T1wb_img, cornersb_img, ipr.iterations, 'Chan-Vese', 'SmoothFactor', ipr.smoothFactor);
            if this.plotmore
                this.plotSegmentation(ac)
                title(sprintf('iterations %i, smooth %g', ipr.iterations, ipr.smoothFactor))
                fn = fullfile(this.anatPath, [this.T1w_ic.fileprefix '_snakes.fig']);
                if ~isfile(fn)
                    savefig(fn)
                end
            end

            % fit back into T1w
            ic = this.T1w_ic.zeros;
            ic.filepath = this.destinationPath;
            ic.fileprefix = [ic.fileprefix '_segmentation'];
            nii = ic.nifti;
            nii.img(this.coords_bb{:}) = ac;
            %%nii.save()
            this.segmentation_ic = mlfourd.ImagingContext2(nii);
        end
        function this = buildCenterlines(this)
            %% builds left and right centerlines, calling this.buildCenterline() for each.
            %  @return this.centerlines_pcs are the pointCloud representation of the centerlines.
            %  @return this.Cs are {L,R} points of the B-spline curve.
            %  @return this.Ps are {L,R} matrices of B-spline control points.

            if ~isempty(this.centerlines_pcs)
                return
            end

            img = logical(this.segmentation_ic);
            imgL = img(1:ceil(this.Nx/2),:,:);
            imgR = zeros(size(img));
            imgR(ceil(this.Nx/2)+1:end,:,:) = img(ceil(this.Nx/2)+1:end,:,:);
            [pcL,CL,PL] = this.buildCenterline(imgL, 'L');
            [pcR,CR,PR] = this.buildCenterline(imgR, 'R');
            this.centerlines_pcs = {pcL pcR};
            this.Cs = {CL CR};
            this.Ps = {PL PR};
        end
        function [pc,C,P] = buildCenterline(this, img, tag)
            %% builds a centerline using mlvg.Hunyadi2021.
            %  @param img is the data upon which a centerline is built.
            %  @param tag is unused.
            %  @return pc is the pointCloud representation of the centerline.
            %  @return C are points of the B-spline curve.
            %  @return P is the matrix of B-spline control points.

            assert(ischar(tag))
            idx = find(img);
            [X,Y,Z] = ind2sub(size(img), idx);             
            M(1,:) = X'; % M are ints cast as double
            M(2,:) = Y';
            M(3,:) = Z';
            this.U = this.NCenterlineSamples; 
            
            P = bspline_estimate(this.k, this.t, M); % double
            C = bspline_deboor(this.k, this.t, P, this.U); % double, ~2x oversampling for Z
            pc = pointCloud(C');
            
            if this.plotmore
                h = figure;
                pcshow(pointCloud(this.T1w_ic, 'thresh', 1000))
                hold on; pcshow(pc.Location, '*m'); hold off;
                fn_fig = fullfile(this.anatPath, sprintf('%s_%s_centerline.fig', this.T1w_ic.fileprefix, tag));
                if ~isfile(fn_fig)
                    saveas(h, fn_fig)
                end
                fn_png = fullfile(this.anatPath, sprintf('%s_%s_centerline.png', this.T1w_ic.fileprefix, tag));
                if ~isfile(fn_png)
                    set(h, 'InvertHardCopy', 'off');
                    set(h,'Color',[0 0 0]); % RGB values [0 0 0] indicates black color
                    saveas(h, fn_png)
                end
            end
            if this.plotdebug
                figure;
                hold all;
                plot3(M(1,:), M(2,:), M(3,:), 'k.');
                plot3(P(1,:), P(2,:), P(3,:), 'b');
                plot3(C(1,:), C(2,:), C(3,:), 'm');
                legend('segmentation', 'control points', 'curve', ...
                    'Location', 'Best');
                hold off;
            end
        end
        function this = buildCORegistrationTargets(this, dyn_ic)
            %% builds CO registration targets comprising time-averaged emissions.
            %  Registration targets are R^3 images.
            
            timeAveraged = dyn_ic.timeAveraged();
            timeAveraged_b25 = timeAveraged.blurred(2.5); % 2.5 mm blurring specified by Fung & Carson
            for i = 1:3
                this.registration.target_ics{i} = timeAveraged_b25; 
            end
        end
        function this = buildRegistrationTargets(this, dyn_ic)
            %% builds registration targets comprising time-averaged emissions for times
            %  sampled at {0.1:0.25,0.25:0.5,0.5:0.75} of maximal whole-brain emissions.
            %  Viz., whole-brain emissions determine the sampling intervals,
            %  but registration targets are R^3 images.

            this.dyn_fileprefix = dyn_ic.fileprefix;
            this.dyn_label = strrep(this.dyn_fileprefix, '_', ' ');
            if contains(dyn_ic.fileprefix, 'CO') || contains(dyn_ic.fileprefix, 'OC')
                this = this.buildCORegistrationTargets(dyn_ic);
                return
            end

            this.wmparc_ic = mlfourd.ImagingContext2(fullfile(this.mriPath, 'wmparc_on_T1w.nii.gz'));
            dyn_avgxyz = dyn_ic.volumeAveraged(logical(this.wmparc_ic));
            dyn_max = dipmax(dyn_avgxyz);
            img = dyn_avgxyz.nifti.img;
            [~,this.it10] = max(img > 0.1*dyn_max);
            [~,this.it25] = max(img > 0.25*dyn_max);
            [~,this.it50] = max(img > 0.5*dyn_max);
            [~,this.it75] = max(img > 0.75*dyn_max);
            
            this.registration.target_ics{1} = dyn_ic.timeAveraged(this.it10:this.it25);
            this.registration.target_ics{2} = dyn_ic.timeAveraged(this.it10:this.it50);
            this.registration.target_ics{3} = dyn_ic.timeAveraged(this.it10:this.it75);
            for i = 1:3
                this.registration.target_ics{i} = this.registration.target_ics{i}.blurred(2.5); % 2.5 mm blurring specified by Fung & Carson
            end
        end
        function [t_idif,ics] = call(this, varargin)
            %% CALL
            %  @param optional toglob, e.g., 'sub-*Dynamic*_on_T1w.nii.gz'
            %  @return table of IDIFs; write table to text file.
            %  @return {ImagingContext2 objects for centerlines}
            
            ip = inputParser;
            addOptional(ip, 'toglob', fullfile(this.petPath, 'sub-*Dynamic*_on_T1w.nii.gz'), @ischar)
            addParameter(ip, 'iterations', 100, @isscalar)
            addParameter(ip, 'smoothFactor', 0, @isscalar)
            parse(ip, varargin{:})
            ipr = ip.Results;

            % build intermediate objects
            this.buildSegmentation(ipr.iterations, 'smoothFactor', ipr.smoothFactor);
            this.buildCenterlines()

            niis = globT(ipr.toglob);
            ics = cell(1, length(niis));
            NIfTI_Filename = cell(1, length(niis));
            tracer = cell(1, length(niis));
            IDIF = cell(1, length(niis));
            alg_ = this.alg;

            for inii = 1:length(niis)

                % sample input function from dynamic PET             
                dyn_ic = mlfourd.ImagingContext2(niis{inii});                
                this.buildRegistrationTargets(dyn_ic)
                this.registerCenterlines('alg', alg_)
                ic = this.pointCloudsToIC();
                ic.filepath = dyn_ic.filepath;
                ic.fileprefix = [dyn_ic.fileprefix '_idifmask'];
                ic.save()
                ics{inii} = ic;
                idif = dyn_ic.volumeAveraged(ic);

                % construct table variables
                NIfTI_Filename{inii} = ic.fqfilename;
                tracer{inii} = this.tracername(ic.fileprefix);
                IDIF{inii} = asrow(this.decay_uncorrected(idif));
            end

            % construct table and write
            t_fileprefix = strsplit(niis{end}, '_on_T1w');
            [~,t_descr] = fileparts(t_fileprefix{1});
            t_idif = table(NIfTI_Filename, tracer, IDIF);
            t_idif.Properties.Description = t_descr;
            t_idif.Properties.VariableUnits = {'', '', 'Bq/mL'};
            t_fqfileprefix = fullfile(this.petPath, [t_descr '_idif']);
            writetable(t_idif, [t_fqfileprefix '.csv']);

            % plot and save
            h = figure;
            hold on
            for irow = 1:length(NIfTI_Filename)
                plot(this.timesMid(tracer{irow}), IDIF{irow}, 'o-')
            end
            xlim([0 350])
            xlabel('time (s)')
            ylabel('activity density (Bq/mL)')
            title('Image-derived Input Functions')
            legend(tracer')
            hold off
            saveas(h, [t_fqfileprefix '.fig'])
            saveas(h, [t_fqfileprefix '.png'])
        end
        function decay_uncorrected = decay_uncorrected(this, idif)
            %  @param idif is an mlfourd.ImagingContext2 containing a double row.
            %  @returns decay_uncorrected, the IDIF as a double row.

            assert(isa(idif, 'mlfourd.ImagingContext2'))
            decay_corrected = idif.nifti.img;
            if contains(idif.fileprefix, 'CO') || contains(idif.fileprefix, 'OC')
                tracer = 'CO';
            end
            if contains(idif.fileprefix, 'Water')
                tracer = 'HO';
            end
            if contains(idif.fileprefix, 'Oxygen')
                tracer = 'OO';
            end
            if contains(idif.fileprefix, 'FDG')
                tracer = 'FDG';
            end
            assert(all(size(decay_corrected) == size(this.taus(tracer))))
            radio = mlpet.Radionuclides(tracer);
            decay_uncorrected = decay_corrected ./ radio.decayCorrectionFactors('taus', this.taus(tracer));
        end
        function [h,h1] = plotRegistered(this, varargin)
            % @param required target, pointCloud.
            % @param required centerlineOnTarget, pointCloud.
            % @param required centerline, pointCloud.
            % @param required laterality, in {'' 'l' 'L' 'r' 'R'}.
            % @return handle(s) for figure(s).
            
            ip = inputParser;
            addRequired(ip, 'target', @(x) isa(x, 'pointCloud'))
            addRequired(ip, 'centerlineOnTarget', @(x) isa(x, 'pointCloud'))
            addRequired(ip, 'centerline', @(x) isa(x, 'pointCloud'))
            addRequired(ip, 'laterality', @(x) ismember(x, {'', 'l', 'L', 'r', 'R'}))
            parse(ip, varargin{:})
            ipr = ip.Results;
            Laterality = upper(ipr.laterality);
        
            h = figure;
            pcshow(ipr.target)
            hold on; 
            pcshow(ipr.centerline.Location, '*g'); 
            pcshow(ipr.centerlineOnTarget.Location, '*m'); 
            hold off;
            title(sprintf('centerline (green -> magenta) on target %s %s', upper(ipr.laterality), this.dyn_label))
            saveas(h, fullfile(this.petPath, sprintf('%s_%s_centerline_target.fig', this.dyn_fileprefix, Laterality)))
            set(h, 'InvertHardCopy', 'off');
            set(h,'Color',[0 0 0]); % RGB values [0 0 0] indicates black color
            saveas(h, fullfile(this.petPath, sprintf('%s_%s_centerline_target.png', this.dyn_fileprefix, Laterality)))
            if this.plotdebug
                h1 = figure;
                pcshowpair(ipr.target, ipr.centerlineOnTarget, 'VerticalAxis', 'Z') % only magenta & green available in R2021b
                title(sprintf('centerline (green) on target (magenta) %s %s', upper(ipr.laterality), this.dyn_label))
                saveas(h1, fullfile(this.petPath, sprintf('%s_%s_centerline_target_mag.fig', this.dyn_fileprefix, Laterality)))
                set(h1, 'InvertHardCopy', 'off');
                set(h1,'Color',[0 0 0]); % RGB values [0 0 0] indicates black color
                saveas(h1, fullfile(this.petPath, sprintf('%s_%s_centerline_target_mag.png', this.dyn_fileprefix, Laterality)))
            end
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
            %% converts point clouds for both hemispheres into ImagingContext objects.
        
            icL = this.pointCloudToIC(this.registration.centerlineOnTarget{1}, varargin{:});
            icL = icL.imdilate(strel('sphere', this.dilationRadius));
            icR = this.pointCloudToIC(this.registration.centerlineOnTarget{2}, varargin{:});
            icR = icR.imdilate(strel('sphere', this.dilationRadius));
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
            this.centerlines_pcs{1} = ...
                this.registerCenterline(this.centerlines_pcs{1}, varargin{:}, 'laterality', 'L');
            this.centerlines_pcs{2} = ...
                this.registerCenterline(this.centerlines_pcs{2}, varargin{:}, 'laterality', 'R');
        end
        function centerlineOnTarget = registerCenterline(this, varargin)
            %  @param required centerline is a pointCloud.
            %  @param optional ic3d is an ImagingContext2 for PET averaged over early times for bolus arrival.
            %         Default is this.registration.target_ics{3}.
            %  @param thresh applies to ic3d.  Default is 25000.
            %  @param alg is from {'ndt', 'icp', 'cpd'}.
            %  @param gridStep preprocesses pcregister* methods.
            %  @param laterality is in {'R' 'L'}.
            %  @return centerlineOnTarget is a pointCloud.
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'centerline', @(x) isa(x, 'pointCloud'))
            addOptional(ip, 'ic3d', this.registration.target_ics{3}, @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'thresh', [], @isnumeric)
            addParameter(ip, 'alg', 'cpd', @(x) ismember(x, {'ndt', 'icp', 'cpd', 'fung'}))
            addParameter(ip, 'gridStep', 1, @isscalar)
            addParameter(ip, 'laterality', '', @(x) ismember(x, {'', 'l', 'L', 'r', 'R'})) % L has indices < Nx/2
            parse(ip, varargin{:})
            ipr = ip.Results;
            ipr.ic3d = this.maskInBoundingBox(ipr.ic3d, ipr.laterality);
            if isempty(ipr.thresh)
                img = ipr.ic3d.nifti.img;
                img = img(img > 0);
                m_ = dipmedian(img);
                s_ = dipstd(img);
                if m_ - s_ > 0
                    ipr.thresh = m_ - s_;
                else
                    ipr.thresh = m_;
                end
            end
            target = pointCloud(ipr.ic3d, 'thresh', ipr.thresh); 
            if this.plotdebug
                figure
                pcshow(target)
                title(sprintf('target %s, thresh %g', upper(ipr.laterality), ipr.thresh))
            end
            centerlineOri = copy(ipr.centerline);
            
            idx = strcmpi(ipr.laterality, 'R') + 1; % idx == 1 <-> left            
            switch ipr.alg
                case 'ndt'
                    [tform,centerlineOnTarget,rmse] = pcregisterndt(centerlineOri, target, ipr.gridStep, ...
                        'Tolerance', [0.01 0.05]);
                case 'icp'
                    if ipr.gridStep ~= 1
                        centerlineOri = pcdownsample(centerlineOri, 'gridAverage', ipr.gridStep);
                    end
                    [tform,centerlineOnTarget,rmse] = pcregistericp(centerlineOri, target, ...
                        'Extrapolate', true, 'Tolerance', [0.01 0.01]);
                case 'cpd'
                    if ipr.gridStep ~= 1
                        centerlineOri = pcdownsample(centerlineOri, 'gridAverage', ipr.gridStep);
                    end
                    [tform,centerlineOnTarget,rmse] = pcregistercpd(centerlineOri, target, ...
                        'Transform', 'Rigid', 'MaxIterations', 100, 'Tolerance', 1e-7); % 'InteractionSigma', 2
                case 'fung'
                    this.registration.tform{idx} = rigid3d(eye(4));
                    rr = mlvg.Reregistration(this.T1w_ic);
                    [tform,centerlineOnTarget,rmse] = rr.pcregistermax( ...
                        this.registration.tform{idx}, centerlineOri, target);
                otherwise
                    error('mlvg:ValueError', ...
                        'Fung2013.registerCenterlines.ipr.alg == %s', ipr.alg)
            end
            this.registration.centerlineOnTarget{idx} = copy(centerlineOnTarget);
            this.registration.tform{idx} = tform;
            this.registration.rmse{idx} = rmse;
       
            if this.ploton
                this.plotRegistered(target, centerlineOnTarget, centerlineOri, ipr.laterality)
            end
        end
        function n = tracername(~, str)
            if contains(str, 'CO')
                n = 'CO';
                return
            end
            if contains(str, 'Oxygen')
                n = 'OO';
                return
            end
            if contains(str, 'Water')
                n = 'HO';
                return
            end
            if contains(str, 'FDG')
                n = 'FDG';
                return
            end
            error('mlvg:ValeError', 'Fung2013.tracername() did not recognize %s', str)
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
        destPath_
        hunyadi_
        projPath_
        subFolder_
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
        function parseDestinationPath(this, dpath)
            assert(contains(dpath, 'CCIR_'), 'Fung2013: destination path must include a project identifier')
            assert(contains(dpath, 'sub-'), 'Fung2013: destination path must include a subject identifier')

            this.destPath_ = dpath;
            ss = strsplit(dpath, filesep);
            [~,idxProjFold] = max(contains(ss, 'CCIR_'));
            this.projPath_ = [filesep fullfile(ss{1:idxProjFold})];
            this.subFolder_ = ss{contains(ss, 'sub-')}; % picks first occurance
        end
    end
end