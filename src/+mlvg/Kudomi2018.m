classdef Kudomi2018 < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
	%% KUDOMI2018  

	%  $Revision$
 	%  was created 13-Jul-2021 23:23:15 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
 	%% It was developed on Matlab 9.10.0.1669831 (R2021a) Update 2 for MACI64.  Copyright 2021 John Joowon Lee.
 	
    properties (Constant)
        ALPHA = log(2)/122.2416
    end
    
	properties
        Acond
        Bcond
        blur
        diagonal_only
        dt
        Nsafeparc % # of voxels in safeparc
        Nskip
        Nx % # of regional samples
        safeparc % mlfourd.ImagingContext2
        scanner
        sessionData
        thresh
        time0
        timeF
    end
    
    properties (Dependent)
        mu
        Nt
        p
        times
        timesMid
        tol
        voxelVolume
    end

	methods
        
        %% GET
        
        function g = get.mu(this)
            g = 1/0.95 + this.ALPHA;
        end
        function g = get.Nt(this)
            g = length(this.times);
        end
        function g = get.p(this)
            %g = 0.8 % Iida 1991, 2000
            g = 1/(1/0.8 + this.ALPHA); % 1/p := 1/lambda + alpha, partition coeff. and radiodecay
        end
        function g = get.times(this)
            g = this.timesMid(1):this.dt:this.timesMid(end);
        end
        function g = get.timesMid(this)            
            timesMid0 = this.scanner.timesMid;
            g = timesMid0(this.time0 <= timesMid0 & timesMid0 <= this.timeF);
        end
        function g = get.tol(this)
            if ~isempty(this.tol_cache)
                g = this.tol_cache;
                return
            end
            g = max(size(this.A))*eps(norm(this.A));
        end
        function g = get.voxelVolume(this)
            g = 1e-3*prod(this.scanner.imagingContext.fourdfp.mmppix); % mL
        end
        
        %%
        
        function arr = aif(this)
            arr = mean(this.aifs(), 1);
        end
        function arr = aifs(this)
            arr = this.X()*this.drho_dt() + this.rho()/this.p;
        end
        function arr = A(this)
            if ~isempty(this.A_cache)
                arr = this.A_cache;
                return
            end
            
            drho_dt_ = this.drho_dt();            
            diag_ = zeros(this.Nx, 1);
            for i = 1:this.Nx
                diag_(i) = sum(drho_dt_(i,:).^2);
            end
            diag_ = (this.Nx - 1)*diag_;            
            offdiag_ = zeros(this.Nx, this.Nx);
            for j = 1:this.Nx
                for i = 1:this.Nx
                    offdiag_(i,j) = -sum(drho_dt_(i,:) .* drho_dt_(j,:));
                end
            end
            
            this.A_cache = diag(diag_) + offdiag_;
            arr = this.A_cache;
            this.Acond = cond(this.A_cache);
        end
        function arr = B(this)
            if ~isempty(this.B_cache)
                arr = this.B_cache;
                return
            end
            
            drho_dt_ = this.drho_dt(); 
            rho_ = this.rho(); 
            diag_ = zeros(this.Nx, 1);
            for i = 1:this.Nx
                diag_(i) = -sum(drho_dt_(i,:).*rho_(i,:));
            end
            diag_ = (this.Nx - 1)*diag_/this.p;
            offdiag_ = zeros(this.Nx, this.Nx);
            for j = 1:this.Nx
                for i = 1:this.Nx
                    offdiag_(i,j) = sum(drho_dt_(i,:) .* rho_(j,:));
                end
            end
            offdiag_ = offdiag_/this.p;
            
            this.B_cache = diag(diag_) + offdiag_;
            arr = this.B_cache;
            this.Bcond = cond(this.B_cache);
        end
        function arr = drho_dt(this)
            if ~isempty(this.drho_dt_cache)
                arr = this.drho_dt_cache;
                return
            end
            smoothed = smoothdata(this.rho(), 'sgolay', 2);
            this.drho_dt_cache = diff(smoothed, 1, 2)/this.dt; % 1st deriv of time
            time_boundary1 = this.drho_dt_cache(:,1);
            this.drho_dt_cache = [time_boundary1 this.drho_dt_cache]; % preserve this.Nt without jitter
            arr = this.drho_dt_cache;
        end
        function arr = d2rho_dt2(this)
            if ~isempty(this.d2rho_dt2_cache)
                arr = this.d2rho_dt2_cache;
                return
            end
            cache = diff(this.rho(), 2, 2)/this.dt^2; % 2nd deriv of time
            time_boundary1 = cache(:,1);
            time_boundaryN = cache(:,end);
            cache = [time_boundary1 cache time_boundaryN]; % preserve this.Nt without jitter
            this.d2rho_dt2_cache = smoothdata(cache, 'movmean', 2);
            arr = this.d2rho_dt2_cache;
        end
        function arr = K1(this)
            arr = 1./diag(this.X());
        end
        function [arr,arr_] = Kvariation(this)
            %% solves \qty[ \mu^2 \rho(x) ]_{N_x \times N_t} = 
            %         \qty[K^{-2}(\vec{x})]_{N_x \times N_x} \qty[ \partial^2_t \rho(x) ]_{N_x \times N_t}. 
            %  @return arr ~ K(\vec{x}).
            %  @return arr_ ~ K^{-2}(\vec{x}).
            
            % solve XA = B
            %arr_ = this.mu^2 * this.rho() .* pinv(this.d2rho_dt2(), this.tol);
            %arr_ = this.mu^2 * this.rho() / this.d2rho_dt2();
            arr_ = lsqminnorm(this.d2rho_dt2()', this.mu^2 * this.rho()')';
            arr = sqrt(1 ./ arr_);
            arr = real(arr);
        end
        function arr = rho(this)
            if ~isempty(this.rho_cache)
                arr = this.rho_cache;
                return
            end
            ic = this.scanner.imagingContext;
            if ~isempty(this.blur)
                ic = ic.blurred(this.blur);
            end
            activityDensityArr = this.contractToNxTacs(ic.fourdfp.img);   
            this.rho_cache = zeros(this.Nx, this.Nt);
            for x = 1:size(this.rho_cache, 1)
                smoothed = smoothdata(activityDensityArr(x,:), 'sgolay');
                interpolated = interp1(this.timesMid, smoothed, this.times, 'makima');
                interpolated(interpolated < 0) = 0;
                interpolated(1:this.Nskip) = 0;
                this.rho_cache(x,:) = interpolated;
            end
            this.rho_cache = this.rho_cache*this.voxelVolume; % Bq
            arr = this.rho_cache;
        end
        function arr = rho_inv(this)
            arr = pinv(this.rho, this.tol);
        end
        function arr = X(this)
            if ~isempty(this.X_cache)
                arr = this.X_cache;
                return
            end
            if this.diagonal_only
                this.X_cache = diag(diag(this.B)./diag(this.A));
            else
                %this.X_cache = this.B .* pinv(this.A, this.tol);
                this.X_cache = this.B/this.A;
                %this.X_cache = lsqminnorm(this.A', this.B')'; % solve XA = B
            end
            arr = this.X_cache;
        end
		  
 		function this = Kudomi2018(varargin)
 			%% KUDOMI2018
            %  @param dt is scalar, in seconds.
            %  @param Nx is scalar, # of voxels to sample.
            %  @param sessionData is an mlpipeline.ISessionData.
 			%  @param tol is numeric.

            ip = inputParser;
            addParameter(ip, 'blur', 9, @isnumeric)
            addParameter(ip, 'diagonal_only', false, @islogical)
            addParameter(ip, 'dt', 2, @isscalar)
            addParameter(ip, 'Nskip', 5, @isscalar)
            addParameter(ip, 'Nx', [], @isnumeric)
            addParameter(ip, 'sessionData', [], @(x) isa(x, 'mlpipeline.ISessionData'))
            addParameter(ip, 'thresh', 0.1, @(x) isscalar(x) && x < 1)
            addParameter(ip, 'time0', [], @isnumeric)
            addParameter(ip, 'timeF', [], @isnumeric)
            addParameter(ip, 'tol', [], @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.blur = ipr.blur;
            this.diagonal_only = ipr.diagonal_only;
            this.dt = ipr.dt;
            this.sessionData = ipr.sessionData;
            this.thresh = ipr.thresh;
 			this.tol_cache = ipr.tol;            
            devkit = mlpet.ScannerKit.createFromSession(this.sessionData);             
            this.scanner = devkit.buildScannerDevice();
            %this.scanner.decayCorrect();
            this.buildSafeparc();
            if ~isempty(ipr.time0)
                this.time0 = ipr.time0;
            else
                this.time0 = this.scanner.timesMid(1);
            end
            if ~isempty(ipr.timeF)
                this.timeF = ipr.timeF;
            else
                this.timeF = this.scanner.timesMid(end);
            end
            this.Nskip = ipr.Nskip;
            if ~isempty(ipr.Nx)
                this.Nx = ipr.Nx;
            else
                this.Nx = this.Nt;
            end
 		end
    end 
    
    %% PROTECTED
    
    properties %(Access = protected)        
 		A_cache % \in \mathcal{R}^2
        B_cache % \in \mathcal{R}^2
        drho_dt_cache % drho(x, t)/dt \in \mathcal{R}^2
        d2rho_dt2_cache % d^2 rho(x, t)/dt^2 \in \mathcal{R}^2
        rho_cache  % drho(x, t)/dt \in \mathcal{R}^2
        tol_cache
        X_cache % \in \mathcal{R}
    end
    
    methods (Access = protected)
        function this = buildSafeparc(this)
            %% builds safeparc, Nsafeparc.
            %  @return this
            
            wmparc1 = this.sessionData.wmparc1OnAtlas('typ', 'mlfourd.ImagingContext2');
            wmparc1_ = wmparc1.fourdfp;
            wmparc1_.img(wmparc1_.img == 1) = 0;
            wmparc1_.img(wmparc1_.img == 6000) = 0;
            this.safeparc = mlfourd.ImagingContext2(wmparc1_);    
            
            safeparc_ = copy(this.safeparc);
            safeparc_ = safeparc_.binarized();
            this.Nsafeparc = safeparc_.dipsum;
        end
        function arr1 = contractToNxTacs(this, img)
            %% bins dynamic img into this.Nx bins ordered by AUCs, masking with this.safeparc and 
            %  thresholding AUCs > 0.1*max(AUCs).
            %  @param img \in \mathcal{R}^{3 + 1}.
            %  @return arr1 \in \mathcal{R}^{1 + 1}; size(arr,1) == this.Nx.
            
            NtF = length(this.timesMid);
            img = img(:,:,:,1:NtF);
            sz = size(img);
            img = reshape(img, [sz(1)*sz(2)*sz(3) sz(4)]); % R^{1 + 1}
            
            % build mask
            bin = this.safeparc.binarized();
            imgmask = bin.fourdfp.img;
            imgmask = reshape(imgmask, [sz(1)*sz(2)*sz(3) 1]); % R^{1 + 1}
            
            % sort aucs, arr
            arr = img(imgmask ~= 0, :);
            aucs = trapz(this.timesMid, arr, 2);
            [aucs_sorted,ind_sorted] = sort(aucs, 'descend');
            arr_sorted = arr(ind_sorted,:);
            
            % apply threshold, keeping higher AUCs
            [~,ind_to_toss] = max(aucs_sorted < this.thresh*max(aucs_sorted));
            if ind_to_toss > 1
                ind_sorted = ind_sorted(1:ind_to_toss);
                arr_sorted = arr_sorted(1:ind_to_toss);
            end
            
            % build binning and average bins
            width = floor(length(ind_sorted)/this.Nx);
            arr1 = zeros(this.Nx, sz(4));
            for i = 1:this.Nx-1
                i0 = width*(i - 1) + 1;
                i1 = width*i;
                arr1(i,:) = mean(arr_sorted(i0:i1,:), 1);
            end
            i0 = width*(this.Nx - 1) + 1;
            arr1(this.Nx,:) = mean(arr_sorted(i0:end,:), 1); % stragglers       
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

