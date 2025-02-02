classdef TimeSeries < handle & mlsystem.IHandle
    %% Explores travelling waves, arousal signal in CO PET
    %  
    %  Created 10-Jul-2024 15:41:00 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 24.1.0.2628055 (R2024a) Update 4 for MACA64.  Copyright 2024 John J. Lee.
    

    properties
        decimate
        hp_thresh_hz  % lower bound ~ 0.01 Hz
        lp_thresh_hz  % higher bound ~ 0.1 Hz

        pet
        t0
    end

    properties (Dependent)
        adj_times
        fileprefix
        fqfp
        frame_duration  % seconds or []
        img
        sampling_freq  % Hz

        iFV  % tentatively just entire 4th ventricle
        brainstem
    end

    methods %% GET
        function g = get.adj_times(this)
            t1 = this.t0 + 1;
            g = this.pet.json_metadata.times(t1:end);            
        end

        function g = get.fileprefix(this)
            g = this.pet.fileprefix;
        end

        function g = get.fqfp(this)
            g = this.pet.fqfp;
        end

        function g = get.frame_duration(this)
            if ~isempty(this.decimate) && isfinite(this.decimate)
                g = this.decimate;
            else
                g = 1;
            end
        end

        function g = get.img(this)
            if ~isempty(this.img_)
                g = this.img_;
                return
            end

            t1 = this.t0 + 1;
            g = this.pet.imagingFormat.img(:,t1:end);
            if ~isempty(this.decimate)
                select = 1:this.decimate:size(g, 2);
                g = g(:, select);
            end
            g = this.scrub1d(g);
            try
                assert(all(isfinite(g), "all"))
                assert(all(g > 0, "all"))
            catch ME
                handexcept(ME)
            end
            this.img_ = g;
        end

        function g = get.iFV(this)
            g = this.img(11, :)';
        end
        function g = get.brainstem(this)
            g = mean(this.img([12, 17, 32], :), 1)';
        end
    end

    methods
        function this = TimeSeries(pet, opts)
            arguments
                pet {mustBeNonempty}
                opts.t0 double = 0
                opts.decimate = []
                opts.hp_thresh_hz = 0.01
                opts.lp_thresh_hz = 0.05
            end

            this.pet = mlfourd.ImagingContext2(pet);
            assert(contains(this.pet.fileprefix, "ParcSchaeffer"))
            this.t0 = opts.t0;
            this.decimate = opts.decimate;
            this.hp_thresh_hz = opts.hp_thresh_hz;
            this.lp_thresh_hz = opts.lp_thresh_hz;
        end

        function dat1 = build_band_passed(this, dat)
            %% Implements butter:  web(fullfile(docroot, 'signal/ref/butter.html?browser=F1help#bucsfmj')) .
            %  See also web(fullfile(docroot, 'signal/ug/practical-introduction-to-digital-filtering.html')) .
            %  Returns:
            %      dat1 same num. type as dat

            if all(dat == 0)
                dat1 = dat;
                return
            end
            if isempty(this.lp_thresh_hz) || isempty(this.hp_thresh_hz) || ...
                    ~isfinite(this.lp_thresh_hz) || ~isfinite(this.hp_thresh_hz)
                dat1 = dat;
                return
            end
            fs = 1/this.frame_duration;
            [z,p,k] = butter(2, [this.hp_thresh_hz, this.lp_thresh_hz]/(fs/2)); % see Matlab doc "Bandpass Butterworth Filter"
            [sos,g] = zp2sos(z, p, k);
            dat1 = filtfilt(sos, g, double(dat));
            if isa(dat, 'single')
                dat1 = single(dat1);
            end
            if isa(dat, 'double')
                dat1 = double(dat1);
            end
        end

        function psi = build_centered(~, psi)
            assert(~isempty(psi))
            if all(psi == 0)
                return
            end
            psi = psi - median(psi, 'all', 'omitnan');
        end

        function psi = build_centered_and_rescaled(this, psi)
            %% Mimics z-score of |psi(t,x)> using median and mad.

            psi = this.build_centered(psi);
            psi = this.build_rescaled(psi);
        end

        function psi = build_global_signal_regressed(this, psi)
            if all(psi == 0)
                return
            end

            psi = psi - median(psi, 2);
        end

        function psi = build_rescaled(~, psi)
            assert(~isempty(psi))
            if all(psi == 0)
                return
            end

            d = mad(abs(psi), 1, 'all');  % flag(mean abs. dev.) ~ 0; flag(median abs. dev.) ~ 1
            psi = psi./d;
        end

        function [img,combo] = parc(this, id, opts)
            %% Provides 7 RSNs from Nick Metcalf and additional parcs derived from Schaeffer or FreeSurfer.
            %  id := "vis" | "sms" | "dan" | "van" | "limbic" | "con" | "dmn" |
            %        "precuneus_grey" | "precuneus_white" |
            %        "grey" | "striatum" | "thalamus" | "white" | "choroid_plexus" | 
            %        "CSF" | "CSF_extracranial" | "ventricles" | 
            %        "iFV" | "brainstem"; case ignored.
            %  opts.laterality := "" |hvr0QER.tza8ncf6tzb containing "l" | containing "r"; case ignored.
            
            %  opts.statistic : is any statistical func taking one argument.
            %
            %  img : Nparc x Ntimes
            %  combo : Nparc x 1

            arguments
                this mlvg.TimeSeries
                id {mustBeText} = ""
                opts.laterality {mustBeTextScalar} = ""
                opts.statistic {mustBeText} = ""
            end

            switch lower(id)
                case "extrinsic"
                    left = 110:182;
                    right = 210:290;
                case "vis"
                    left = 110:123;
                    right = 210:224;
                case "sms"
                    left = 124:139;
                    right = 225:243;
                case "dan"
                    left = 140:152;
                    right = 244:256;
                case "van"
                    left = 153:163;
                    right = 257:267;
                case "limbic"
                    left = 164:169;
                    right = 268:273;
                case "con"
                    left = 170:182;
                    right = 274:290;
                case "dmn"
                    left = 183:209;
                    right = 291:309;
                case "precuneus_grey"
                    left = 205:208;
                    right = 307:309;
                case "precuneus_white"
                    left = 65;
                    right = 99;
                case "grey"
                    left = 110:209;
                    right = 210:309;
                case "striatum"  % caudate and lentiform nuclei
                    left = 7:9;
                    right = 26:28;
                case "thalamus"
                    left = 6;
                    right = 25;
                case "hippocampus"
                    left = 13;
                    right = 29;
                case "amygdala"
                    left = 14;
                    right = 30;
                case "white"
                    left = 42:75;
                    right = 76:109;
                case "choroid_plexus" 
                    left = 19;
                    right = 33;
                case "csf"  % ?
                    left = 15;
                    right = left;  % provides correct result for inadvertent laterality requests of unary parcs
                case "csf_extracranial" 
                    left = 40;
                    right = left;
                case "ventricles"  % L lat, L inf lat, 3rd, 4th, R lat, R inf lat
                    left = [2, 3, 10, 11, 21, 22];
                    right = left;
                case "arousal"
                    left = [11, 12];
                    right = left;
                case "ifv" 
                    left = 11;
                    right = left;
                case "brainstem"  % including ventral diencephalon
                    left = [12, 17, 32];
                    right = left;
                case "cerebellum"
                    left = [4, 5];
                    right = [23, 24];
                case "cerebellar_cortex"
                    left = 5;
                    right = 24;
                case "cerebellar_white"
                    left = 4;
                    right = 23;
                case "vermis"
                    left = 38;
                    right = left;
                case "skull"
                    left = 37;
                    right = left;
                case "air"
                    left = 36;
                    right = left;
                otherwise
                    error("mlvg:ValueError", stackstr()+":lower(id)->"+lower(id))
            end

            % assemble elemental parcs
            img = this.img;
            if contains(opts.laterality, "l",  IgnoreCase=true) 
                combo = left;
                img = img(left, :);
            elseif contains(opts.laterality, "r",  IgnoreCase=true)
                combo = right;
                img = img(right, :);
            else
                if length(left) ~= length(right) || ~all(left == right)
                    combo = [left, right];
                    img = img([left, right], :);
                else  % unary parc
                    combo = left;
                    img = img(left, :);
                end
            end

            % apply statistic
            if ~isemptytext(opts.statistic)
                img = eval(opts.statistic+"(img, 1)");
            end
        end

        function h = plot(this, img, opts)
            %% 
            %  img double
            %  opts.id {mustBeText} = ""
            %  opts.use_cbrewer2 logical = true
            %  opts.colormap {mustBeText} = "Accent"
            %
            %  For cbrewer2, see also 
            %  https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/58350/versions/2/screenshot.png

            arguments
                this mlvg.TimeSeries
                img double
                opts.id {mustBeText} = ""
                opts.combo {mustBeNumeric} = []
                opts.use_cbrewer2 logical = true
                opts.colormap {mustBeText} = "Accent"
                opts.new_figure logical = false
                opts.LineWidth double = 0.5
                opts.title {mustBeText} = ""
            end
            if ~isemptytext(opts.id)
                opts.id = opts.id + " ";
            end
            if isemptytext(opts.title)
                opts.title = this.pet.fileprefix+" : "+opts.id;
            end

            % color mapping
            Nparc = size(img, 1);
            if opts.use_cbrewer2
                cmap = eval(sprintf("cbrewer2('%s', %i)", char(opts.colormap), Nparc));
            else
                cmap = eval(sprintf("%s(%i)", opts.colormap, Nparc));
            end
            assert(isnumeric(cmap))

            % plot to figure; add colorbar     
            if opts.new_figure
                h = figure;
            else
                h = [];
            end
            hold on;
            for p = 1:Nparc
                plot(img(p, :), Color=cmap(p, :), LineWidth=opts.LineWidth)
            end
            hold off;
            colormap(cmap);
            if Nparc >= 9
                c = colorbar(Direction="reverse");
                c.Label.String = opts.id+"Schaeffer parcel";
                clim([1, Nparc]);                
            end

            % annotate
            ylabel("activity (Bq/mL)"); 
            xlabel("time of mid-frame (sec)");
            if Nparc < 9 && ~isempty(opts.combo)
                legend("parc = "+ascol(opts.combo))
            end
            % fontsize(scale=2);
            title(opts.title, FontSize=9, Interpreter="none")
        end

        function s = size(this, varargin)
            s = size(this.img, varargin{:});
        end
    end

    methods (Static)
        function img = scrub1d(img)
            %% img ~ parcs x times

            Nt = size(img, 2);
            defects = ~isfinite(img) | ~(img > 0);
            [p,t] = find(defects);
            assert(length(p) == length(t))
            for idx = 1:length(p)
                if t(idx) == 1
                    img(p(idx), 1) = img(p(idx), 2);
                    continue
                end
                if t(idx) == Nt
                    img(p(idx), Nt) = img(p(idx), Nt - 1);
                    continue
                end
                img(p(idx), t(idx)) = mean([img(p(idx), t(idx) - 1), img(p(idx), t(idx) + 1)]);
            end
        end
        function result = scrub2d(matrix)
            %% https://claude.ai/chat/c053dfc0-e6e1-41f3-8931-784c93501f10

            % Create a kernel for neighboring elements (including diagonals)
            kernel = ones(3,3);
            kernel(2,2) = 0;  % Exclude the center element

            % Count non-zero neighbors
            nonZeroCount = conv2(double(matrix ~= 0), kernel, 'same');

            % Sum of neighboring elements
            neighborSum = conv2(matrix, kernel, 'same');

            % Calculate average of non-zero neighbors
            avgNeighbors = neighborSum ./ nonZeroCount;

            % Replace zeros with average of non-zero neighbors
            result = matrix;
            zeroMask = (matrix == 0);
            result(zeroMask) = avgNeighbors(zeroMask);

            % Handle edge cases where all neighbors are zero
            allZeroNeighbors = (nonZeroCount == 0);
            result(allZeroNeighbors) = 0;
        end
    end

    %% PRIVATE

    properties (Access=private)
        img_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
