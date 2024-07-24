classdef TimeSeries < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 10-Jul-2024 15:41:00 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 24.1.0.2628055 (R2024a) Update 4 for MACA64.  Copyright 2024 John J. Lee.
    

    properties
        pet
        t0
    end

    properties (Dependent)
        adj_times
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

        function g = get.img(this)
            if ~isempty(this.img_)
                g = this.img_;
                return
            end

            t1 = this.t0 + 1;
            this.img_ = this.pet.imagingFormat.img(:,t1:end);
            g = this.img_;
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
            end

            this.pet = mlfourd.ImagingContext2(pet);
            assert(contains(this.pet.fileprefix, "ParcSchaeffer"))
            this.t0 = opts.t0;
        end

        function [img,combo] = parc(this, id, opts)
            %% Provides 7 RSNs from Nick Metcalf and additional parcs derived from Schaeffer or FreeSurfer.
            %  id := "vis" | "sms" | "dan" | "van" | "limbic" | "con" | "dmn" |
            %        "precuneus_grey" | "precuneus_white" |
            %        "grey" | "striatum" | "thalamus" | "white" | "choroid_plexus" | 
            %        "CSF" | "CSF_extracranial" | "ventricles" | 
            %        "iFV" | "brainstem"; case ignored.
            %  opts.laterality := "" | containing "l" | containing "r"; case ignored.
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
                    left = 204:208;
                    right = 307:309;
                case "precuneus_white"
                    left = 65;
                    right = 99;
                case "grey"
                    left = 110:209;
                    right = 210:309;
                case "striatum"  % caudate and lentiform nuclei
                    left = 11:13;
                    right = 50:52;
                case "thalamus"
                    left = 6;
                    right = 25;
                case "white"
                    left = 42:75;
                    right = 76:109;
                case "choroid_plexus" 
                    left = 19;
                    right = 33;
                case "CSF"  % ?
                    left = 15;
                    right = left;  % provides correct result for inadvertent laterality requests of unary parcs
                case "CSF_extracranial" 
                    left = 40;
                    right = left;
                case "ventricles"  % L lat, L inf lat, 3rd, 4th, R lat, R inf lat
                    left = [2, 3, 10, 11, 21, 22];
                    right = left;
                case "iFV" 
                    left = 11;
                    right = left;
                case "brainstem"  % including ventral diencephalon
                    left = [12, 17, 32];
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
    end

    %% PRIVATE

    properties (Access=private)
        img_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
