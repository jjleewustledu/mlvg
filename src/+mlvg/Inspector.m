classdef Inspector < handle
    %% line1
    %  line2
    %  
    %  Created 12-Aug-2025 16:08:43 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.
    
    properties (Constant)
        WATER_DENSITY = 0.9982 % pure water at 20 C := 0.9982 mL/g; tap := 0.99823
        BLOOD_DENSITY = 1.06 % https://hypertextbook.com/facts/2004/MichaelShmukler.shtml; human whole blood 37 C
        BRAIN_DENSITY = 1.05 % Torack et al., 1976
        LC = 0.81
        PLASMA_DENSITY = 1.03
    end

    properties
        N_parcels = 309
        N_schaefer = 200
    end

    properties (Dependent)
        indices_wb
        indices_gm
        indices_wm  % no cerebellum
        indices_brainstem
        indices_cerebellum
        indices_csf
        indices_subcortex  % no cerebellum
        indices_schaef  % 309 non-zero elements
        N_symmetric  % averaging LH with RH, e.g., 209
        N_notcort  % not cortical, e.g., 109
    end

    methods  %% GET/SET
        function g = get.indices_wb(this)
            g = [this.indices_gm, this.indices_wm, this.indices_subcortex, this.indices_cerebellum, this.indices_brainstem];
        end
        function g = get.indices_gm(this)
            g = 20001:20200;
        end
        function g = get.indices_wm(this)
            g = [2, 41, 3001:3003, 3005:3035, 4001:4003, 4005:4035];
        end
        function g = get.indices_subcortex(this)
            g = [ ...
                10, 11, 12, 13, 17, 18, ...
                26, 28, 43, 44, 49, 50, 51, 52, ...
                53, 54, 58, 60, 174];
        end
        function g = get.indices_brainstem(this)
            g = 16;
        end
        function g = get.indices_cerebellum(this)
            g = [7, 8, 46, 47, 172];
        end
        function g = get.indices_csf(this)
            g = [4, 5, 14, 15, 24, 31, 43, 44, 63, 257];
        end
        function g = get.indices_schaef(this)
            if ~isempty(this.indices_schaef_)
                g = this.indices_schaef_;
                return
            end

            ld = load(fullfile(getenv("SINGULARITY_HOME"), "CCIR_01211", "indices_schaef.mat"));
            g = ld.indices_schaef;
            g = g(g ~= 0);  % assurance
            this.indices_schaef_ = g;
        end
        function g = get.N_symmetric(this)
            g = this.N_parcels - this.N_schaefer/2;
        end
        function g = get.N_notcort(this)
            g = this.N_parcels - this.N_schaefer;
        end
    end

    methods
        function this = Inspector(varargin)
            warning("off", "MATLAB:table:ModifiedAndSavedVarnames")
        end

        function inspect_aging_cmrglc(this, fdg)
            arguments
                this mlvg.Inspector
                fdg {mustBeText}  % array of fqfn
            end

            measure_name = "CMR_{glc}";
            measure_units = "\mu mol hg^{-1}min^{-1}";
            converter = @(x) x * 6000 / this.BRAIN_DENSITY;

            data_k1_k2 = this.load(fdg, measure_index=1);
            data_k2 = this.load(fdg, measure_index=2);
            data_k1 = data_k1_k2 .* data_k2;
            data_k3_k4 = this.load(fdg, measure_index=3);
            data_k4 = this.load(fdg, measure_index=4);
            data_k3 = data_k3_k4 .* data_k4;
            data_v1 = this.load_matching_v1(fdg);
            assert(length(data_v1) == length(data_k1))

            denom = data_k2 + data_k3;
            data_Ki = data_v1 .* data_k1 .* data_k3 ./ denom;
            glc = this.load_matching_glc(fdg);
            data = data_Ki .* asrow(glc) / this.LC;
            data = converter(data);
            dstruct = this.weighted_average(data, fdg);  % 309 x Nses =: struct with fields wb, dmn, gm, wm, subcortex, each ~ 1 x Nses
            age = this.load_matching_age(fdg);  % ~ 1 x Nses

            %% plotting trends

            this.plot_age_vs_measurements_combined(age, dstruct, ...
                title="", ...
                ylabel=sprintf("%s (%s)", measure_name, measure_units));
        end

        function inspect_aging_cmro2(this, oo)
            arguments
                this mlvg.Inspector
                oo {mustBeText}  % array of fqfn
            end

            measure_name = "CMRO_2";
            measure_units = "\mu mol hg^{-1}min^{-1}";
            converter = @(x) x * 6000 / this.BRAIN_DENSITY;

            data_oef = this.load(oo, measure_index=1);
            data_v_post = this.load(oo, measure_index=3);
            data_f = this.load_matching_f(oo);
            assert(length(data_f) == length(data_oef))

            o2_content = this.load_matching_o2_content(oo);
            data = data_f .* data_oef .* (data_v_post / 0.835) .* o2_content;
            data = converter(data);
            dstruct = this.weighted_average(data, oo);  % 309 x Nses =: struct with fields wb, dmn, gm, wm, subcortex, each ~ 1 x Nses
            age = this.load_matching_age(oo);  % ~ 1 x Nses

            %% plotting trends

            this.plot_age_vs_measurements_combined(age, dstruct, ...
                title="", ...
                ylabel=sprintf("%s (%s)", measure_name, measure_units));
        end


        %% UTILITIES

        function fig = plot_age_vs_measurements_combined(this, a, dstruct, viridis_func, opts)
            % PLOT_AGE_VS_MEASUREMENTS_COMBINED Plot all measurements vs age with trend lines
            %
            % Inputs:
            %   a - Age vector (length L)
            %   w, x, y, z - Measurement vectors (same length as a)
            %   viridis_func - Optional function handle for viridis colormap
            %
            % Output:
            %   fig - Figure handle
            %
            % Example:
            %   fig = plot_age_vs_measurements_combined(age, w_data, x_data, y_data, z_data, @viridis);

            arguments
                this mlvg.Inspector
                a {mustBeNumeric}
                dstruct struct
                viridis_func function_handle = @viridis
                opts.title {mustBeTextScalar} = ""
                opts.ylabel {mustBeTextScalar} = "Measurements"
            end
            u = dstruct.cerebellum;
            v = dstruct.brainstem;
            w = dstruct.wb;
            x = dstruct.gm;
            y = dstruct.wm;
            z = dstruct.subcortex;

            % Create figure
            fig = figure('Position', [100, 100, 1200, 800]);
            hold on;

            % Get viridis colors
            cmap = viridis_func();

            % Select 6 distinct colors from viridis colormap
            n_colors = size(cmap, 1);
            color_indices = round(linspace(0.2*n_colors, 0.9*n_colors, 6));
            colors = cmap(color_indices, :);

            % Data vectors and labels
            data_vectors = {w, x, y, z, u, v};
            labels = {'wb', 'gm', 'wm', 'subcortex', 'cerebellum', 'brainstem'};
            markers = {'^', 'o', 'o', 'o', 'o', 'o'};  % triangle for wb, circles for others

            % Point size and transparency
            point_size = 100;
            alpha = 0.5;

            % Store handles for legend
            scatter_handles = [];
            trend_handles = [];
            legend_labels = {};

            % Process each data vector
            for idx = 1:6
                data = data_vectors{idx};
                label = labels{idx};
                color = colors(idx, :);
                marker = markers{idx};

                % Convert to column vectors
                age_vec = a(:);
                data_vec = data(:);

                % Remove NaN values
                valid_mask = ~(isnan(age_vec) | isnan(data_vec));
                age_clean = age_vec(valid_mask);
                data_clean = data_vec(valid_mask);

                % Scatter plot with transparency
                h_scatter = scatter(age_clean, data_clean, point_size, color, marker, 'filled');
                h_scatter.MarkerFaceAlpha = alpha;
                h_scatter.MarkerEdgeColor = 'none';

                % Store handle and create label
                scatter_handles(end+1) = h_scatter;
                legend_labels{end+1} = sprintf('%s (n=%d)', label, length(age_clean));

                % Calculate and plot trend line if we have enough data
                if length(age_clean) > 1
                    % Linear regression using polyfit
                    p = polyfit(age_clean, data_clean, 1);
                    slope = p(1);
                    intercept = p(2);

                    % Calculate R-squared
                    y_fit = polyval(p, age_clean);
                    ss_res = sum((data_clean - y_fit).^2);
                    ss_tot = sum((data_clean - mean(data_clean)).^2);
                    r_squared = 1 - (ss_res / ss_tot);

                    % Generate trend line
                    x_trend = linspace(min(age_clean), max(age_clean), 100);
                    y_trend = polyval(p, x_trend);

                    % Plot trend line
                    h_trend = plot(x_trend, y_trend, '--', ...
                        'Color', color, ...
                        'LineWidth', 1.5);

                    % Store handle and create label
                    trend_handles(end+1) = h_trend;
                    legend_labels{end+1} = sprintf('%s trend (R²=%.3f)', label, r_squared);
                end
            end
    
            % Calculate data range for y-axis limits with padding
            all_data = [w(:); x(:); y(:); z(:); u(:); v(:)];
            all_data_clean = all_data(~isnan(all_data));

            if ~isempty(all_data_clean)
                data_min = min(all_data_clean);
                data_max = max(all_data_clean);
                data_range = data_max - data_min;

                % Add 10% padding to top and bottom
                padding = 0.15 * data_range;
                ylim([data_min - padding, data_max + padding]);
            end

            % Calculate x-axis limits with padding
            all_ages = a(:);
            all_ages_clean = all_ages(~isnan(all_ages));

            if ~isempty(all_ages_clean)
                age_min = min(all_ages_clean);
                age_max = max(all_ages_clean);
                age_range = age_max - age_min;

                % Add 5% padding to left and right
                x_padding = 0.05 * age_range;
                xlim([age_min - x_padding, age_max + x_padding]);
            end

            % Formatting
            xlabel('Age (years)', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel(opts.ylabel, 'FontSize', 12, 'FontWeight', 'bold');
            title(opts.title, 'FontSize', 14, 'FontWeight', 'bold');

            % Grid
            grid on;
            grid minor;
            set(gca, 'GridAlpha', 0.3);
            set(gca, 'MinorGridAlpha', 0.15);

            % Legend
            all_handles = [scatter_handles, trend_handles];
            legend_labels = legend_labels([1, 3, 5, 7, 9, 11, 2, 4, 6, 8, 10, 12]);
            legend(all_handles, legend_labels, ...
                'Location', 'best', ...
                'FontSize', 10, ...
                'NumColumns', 2);

            % Improve axes appearance
            set(gca, 'FontSize', 10);
            box on;

            hold off;

            fontsize(scale=1.6)
        end

        function cmap = viridis_colormap(this)
            % VIRIDIS_COLORMAP Generate a viridis-like colormap
            %
            % This creates a colormap similar to Python's viridis
            % Returns a 256x3 matrix of RGB values

            % Define key colors from viridis palette
            key_colors = [
                0.267, 0.005, 0.329;  % Dark purple
                0.283, 0.141, 0.458;  % Purple
                0.254, 0.265, 0.530;  % Blue-purple
                0.207, 0.371, 0.553;  % Blue
                0.164, 0.471, 0.558;  % Blue-green
                0.128, 0.567, 0.551;  % Green-blue
                0.135, 0.659, 0.518;  % Green
                0.267, 0.749, 0.441;  % Yellow-green
                0.478, 0.821, 0.318;  % Yellow-green
                0.741, 0.873, 0.150;  % Yellow
                0.993, 0.906, 0.144;  % Bright yellow
                ];

            % Interpolate to create 256 colors
            n_colors = 256;
            x = linspace(1, size(key_colors, 1), size(key_colors, 1));
            xi = linspace(1, size(key_colors, 1), n_colors);

            cmap = zeros(n_colors, 3);
            for channel = 1:3
                cmap(:, channel) = interp1(x, key_colors(:, channel), xi, 'pchip');
            end

            % Ensure values are in [0, 1] range
            cmap = max(0, min(1, cmap));
        end

        function heatmap(this, data, opts)
            arguments
                this mlvg.Inspector
                data {mustBeNumeric}
                opts.measure {mustBeText} = "cmro2"
                opts.measure_name {mustBeTextScalar} = "CMRO2";
                opts.measure_units {mustBeTextScalar} = "\mu mol min^{-1} hg{-1}";
            end

            % Sort sensors by median or mean to highlight the clustering
            medians = median(data, 2);
            [~, sortIdx] = sort(medians, 'descend');
            data_sorted = data(sortIdx, :);

            % Create heatmap
            figure('Position', [100, 100, 800, 1200]);
            imagesc(data_sorted');
            colorbar;
            xlabel('Schaefer parcel (sorted)');
            ylabel(sprintf('%s (%s)', opts.measure_name, opts.measure_units));
            title(sprintf('%s Heatmap', opts.measure_name));
        end

        function raincloud(this, data, opts)
            arguments
                this mlvg.Inspector
                data {mustBeNumeric}
                opts.measure {mustBeText} = "cmro2"
                opts.measure_name {mustBeTextScalar} = "CMRO2";
                opts.measure_units {mustBeTextScalar} = "\mu mol min^{-1} hg{-1}";
                opts.xlabel {mustBeTextScalar} = "[^{15}O] PET session"
            end

            M = size(data, 1);
            figure('Position', [100, 100, 1200, 1400]);

            % Plot quantiles as thin lines
            subplot(1, 2, 1);
            quantiles = prctile(data, [10 25 50 75 90], 2);
            hold on;
            for m = 1:M
                color = this.parc_color(m, M);

                % Plot median as thicker line
                plot(quantiles(m, 3), m, 'o', 'Color', color, 'MarkerSize', 2);
                % Plot IQR
                plot([quantiles(m, 2), quantiles(m, 4)], [m, m], '-', ...
                    'Color', color, 'LineWidth', 1.5);
                % Plot 10-90 percentile
                plot([quantiles(m, 1), quantiles(m, 5)], [m, m], '-', ...
                    'Color', color, 'LineWidth', 0.5);
            end
            if strcmpi(opts.measure, "k_3")
                xlim([0, 20 / 60])
            end
            if strcmpi(opts.measure, "k_3/k_4")
                xlim([0, 50 / 60])
            end
            ylim([0.5, M+0.5]);  % Add 0.5 padding to match imagesc behavior
            set(gca, 'YDir', 'reverse');
            ylabel('Schaefer parcel');
            xlabel(sprintf('%s (%s)', opts.measure_name, opts.measure_units));
            title(sprintf('%s quantiles', opts.measure_name));

            % Add a compact heatmap showing density
            subplot(1, 2, 2);
            imagesc(data);
            colormap(gca, viridis);  % Apply viridis only to current axes
            ylabel('Schaefer parcel');
            xlabel(opts.xlabel)
            c = colorbar;
            clim([dipmin(data), dipmax(data)]);
            ylabel(c, sprintf('%s (%s)', opts.measure_name, opts.measure_units), 'Rotation', 270);
            title('Raw Data Heatmap');

            fontsize(scale=1.6)
        end

        function [data,found] = load(this, fqfns, opts)
            %% loads fqfns into data matrix, ignoring fqfns not on filesystem.  
            %  Returned data selects measurement index, applies numerical conversions, and
            %  averages cortical lhs with rhs.  
            %  Returns:
            %      data ~ 209 x N_ses

            arguments
                this mlvg.Inspector
                fqfns {mustBeText}
                opts.measure_index {mustBeScalarOrEmpty} = 1  % idx of parameter in Dynesty/Idif2025 results
                opts.converter function_handle = @(x) x  % rescaling units, etc.
            end

            %fqfns = fqfns(isfile(fqfns));

            found = asrow(true(size(fqfns)));
            measures = [];
            fidx = 0;
            for fqfn = fqfns
                fidx = fidx + 1;
                if ~isfile(fqfn)
                    found(fidx) = false;
                    measures = [measures, nan(this.N_parcels, 1)]; %#ok<AGROW>
                    continue
                end
                ifc = mlfourd.ImagingFormatContext2(fqfn);
                N_parcels_ = size(ifc.img, 1);
                if N_parcels_ ~= numel(this.indices_schaef)
                    fprintf("%s: %s has %g parcels\n", stackstr(), ifc.filename, N_parcels_);
                    ifc = this.repair_schaefer(ifc);
                end
                measure = ascol(opts.converter(ifc.img(:, opts.measure_index)));
                measures = [measures, measure];  %#ok<AGROW> % N_parcels x N_ses
            end

            data = this.symmetrize_cortex(measures);
        end

        function data = load_matching_age(~, fqfns)
            %% Returns:
            %      data ~ 1 x N_ses

            % obj supplying glc conc ~ $\mu$mol/mL
            clin = mlvg.Clinical();

            data = nan(size(fqfns));
            for fidx = 1:length(fqfns)
                s = mlpipeline.Bids.filename2struct(fqfns(fidx));
                datestr = extractAfter(s.ses, "ses-");
                pet_date = datetime(datestr, InputFormat="yyyyMMddHHmmss");
                try
                    data(fidx) = clin.sub_to_age(s.sub, pet_date);
                catch
                    data(fidx) = nan;
                end
            end
        end

        function data = load_matching_f(this, fqfns)
            %% Returns:
            %      data ~ 209 x N_ses, mL/cm^3

            fqfns_f = [];
            for fqfn = fqfns
                pth = fileparts(fqfn);
                pth_sub = extractBefore(pth, filesep + "ses-");
                % sub-108007_ses-20210219145054_trc-ho_proc-ParcSchaeffer-reshape-to-schaeffer-schaeffer-finite-TissueIO-Boxcar--Raichle1983-qm.nii.gz
                glob = mglob(fullfile(pth_sub, "ses-*", "pet", "sub-*_ses-*_trc-ho_proc-ParcSchaeffer*Boxcar--Raichle1983-qm.nii.gz"));
                if isempty(glob)
                    fqfn_f = "";
                else
                    fqfn_f = glob(1);
                end
                fqfns_f = [fqfns_f, fqfn_f]; %#ok<AGROW>
            end

            data = this.load(fqfns_f);
        end

        function data = load_matching_glc(~, fqfns)
            %% Returns:
            %      data ~ 1 x N_ses

            % obj supplying glc conc ~ $\mu$mol/mL
            clin = mlvg.Clinical();

            data = nan(size(fqfns));
            for fidx = 1:length(fqfns)
                s = mlpipeline.Bids.filename2struct(fqfns(fidx));
                try
                    data(fidx) = clin.sub_to_glc(s.sub);
                catch
                    data(fidx) = nan;
                end
            end
        end

        function data = load_matching_o2_content(~, fqfns)
            %% Returns:
            %      data ~ 1 x N_ses

            % obj supplying glc conc ~ $\mu$mol/mL
            clin = mlvg.Clinical();

            data = nan(size(fqfns));
            for fidx = 1:length(fqfns)
                s = mlpipeline.Bids.filename2struct(fqfns(fidx));
                try
                    data(fidx) = clin.sub_to_o2_content(s.sub);
                catch
                    data(fidx) = nan;
                end
            end
        end

        function data = load_matching_v1(this, fqfns)
            %% Returns:
            %      data ~ 209 x N_ses, mL/cm^3/s

            fqfns_v1 = [];
            for fqfn = fqfns
                pth = fileparts(fqfn);
                pth_sub = extractBefore(pth, filesep + "ses-");
                glob = mglob(fullfile(pth_sub, "ses-*", "pet", "sub-*_ses-*_martinv1.nii.gz"));
                if isempty(glob)
                    fqfn_v1 = "";
                else
                    fqfn_v1 = glob(1);
                end
                fqfns_v1 = [fqfns_v1, fqfn_v1]; %#ok<AGROW>
            end

            data = this.load(fqfns_v1);
        end        

        function color = parc_color(this, parc_idx, N_idx)
            switch N_idx
                case {209,309}
                    if parc_idx > this.N_notcort
                        color = [0.5, 0.5, 0.7];  % blue-grey for cortex
                    elseif 42 <= parc_idx && parc_idx <= this.N_notcort  % light grey for superficial white
                        color = [0.7, 0.7, 0.7];
                    elseif parc_idx == 1 || parc_idx == 20  % medium grey for deep white
                        color = [0.5, 0.5, 0.5];
                    else
                        color = [0.8, 0.3, 0.2];  % maroon for subcortex & white
                    end
                case {210,310}
                    if parc_idx > this.N_notcort + 1
                        color = [0.5, 0.5, 0.7];  % blue-grey for cortex
                    elseif 43 <= parc_idx && parc_idx <= this.N_notcort + 1  % light grey for superficial white
                        color = [0.7, 0.7, 0.7];
                    elseif parc_idx == 36
                        color = [1, 0, 0];  % 80 non-WM-hypointensities 164 108 226 0 2
                    elseif parc_idx == 1 || parc_idx == 20  % medium grey for deep white
                        color = [0.5, 0.5, 0.5];
                    else
                        color = [0.8, 0.3, 0.2];  % maroon for subcortex & white
                    end
                otherwise
                    error("mlvg:ValueError", stackstr());
            end
        end
    
        function ifc = registered_schaefer(~, fqfn)
            % sub-108121_ses-20231030144329_trc-fdg_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames-schaeffer.nii.gz
            g = mglob(fullfile(fileparts(fqfn), "sub-*_ses-*BrainMoCo2-createNiftiMovingAvgFrames-schaeffer.nii.gz"));
            assert(~isempty(g))
            ifc = mlfourd.ImagingFormatContext2(g(1));
        end

        function ifc = repair_schaefer(this, ifc, opts)
            arguments
                this mlvg.Inspector
                ifc mlfourd.ImagingFormatContext2
                opts.do_replace_file logical = false
            end

            try
                registered_schaef = this.registered_schaefer(ifc.fqfn);  % sometimes fails
                indices_registered = unique(registered_schaef.img);
                indices_registered = setdiff(indices_registered, 0);  % exclude 0
                assert(~isempty(indices_registered))

                % in case of superfluous call to repair
                if numel(this.indices_schaef) == numel(indices_registered)
                    return
                end

                img = ifc.img;
                nans_1xwidth = nan(1, size(img, 2));

                % excess indices_registered
                if numel(this.indices_schaef) < numel(indices_registered)
                    [excess,poss_e] = setdiff(indices_registered, this.indices_schaef);
                    for pos_e = sort(poss_e)
                        img(pos_e,:) = [];
                    end
                end

                % missing indices_registered
                if numel(indices_registered) < numel(this.indices_schaef)
                    [missing,poss_m] = setdiff(this.indices_schaef, indices_registered);
                    for pos_m = poss_m
                        img = ...
                            [img(1:pos_m-1, :); ...
                            nans_1xwidth; ...
                            img(pos_m:end, :)];
                    end
                end

                ifc.img = img;
                if opts.do_replace_file
                    ifc.save();
                end
            catch ME
                ifc.img = nan(309, size(ifc.img, 2));
                handwarning(ME.message)
                return
            end
        end

        function mat = repair_schaefer_mat(this, mat, indices_registered)
            arguments
                this mlvg.Inspector
                mat {mustBeNumeric}
                indices_registered {mustBeNumeric}
            end
            assert(all(size(indices_registered) == size(mat)))

            try
                % in case of superfluous call to repair
                if numel(this.indices_schaef) == numel(indices_registered)
                    return
                end

                nans_1xwidth = nan(1, size(mat, 2));

                % excess indices_registered
                if numel(this.indices_schaef) < numel(indices_registered)
                    [excess,poss_e] = setdiff(indices_registered, this.indices_schaef);
                    for pos_e = sort(poss_e)
                        mat(pos_e,:) = [];
                    end
                end

                % missing indices_registered
                if numel(indices_registered) < numel(this.indices_schaef)
                    [missing,poss_m] = setdiff(this.indices_schaef, indices_registered);
                    for pos_m = poss_m
                        mat = ...
                            [mat(1:pos_m-1, :); ...
                            nans_1xwidth; ...
                            mat(pos_m:end, :)];
                    end
                end
            catch ME
                mat = nan(309, size(mat, 2));
                handwarning(ME.message)
                return
            end
        end
        
        function data = symmetrize_cortex(this, measures)
            %% 309 x Nses =: 209 x Nses

            % abbrev.
            Nnc = this.N_notcort;
            Nsplit = Nnc + this.N_schaefer / 2;

            data = nan(this.N_symmetric, size(measures, 2));  % average cortical lhs, rhs
            data(Nnc+1:end, :) = (measures(Nnc+1:Nsplit, :) + measures(Nsplit+1:this.N_parcels, :)) / 2;
            data(1:Nnc, :) = measures(1:Nnc, :);  % wm and subcortical
        end

        function w = weights_all_schaefer(this, fqfn)
            %% ensures w ~ 309 x 1, \sum_a w^a = 1.

            registered_schaef = this.registered_schaefer(fqfn);
            indices_registered = unique(registered_schaef.img);
            indices_registered = setdiff(indices_registered, 0);  % exclude 0
            
            w = nan(size(indices_registered));  % column
            for idx = 1:length(indices_registered)
                w(idx) = sum(indices_registered(idx) == registered_schaef.img, "all");
            end
            w = w / sum(w, "all");
            if length(w) ~= length(this.indices_schaef)
                w = this.repair_schaefer_mat(w, indices_registered);
            end
        end

        function dstruct = weighted_average(this, data, pet_fqfn)
            %% 309 x Nses =: struct with fields wb, gm, dmn, wm, subcortex, each ~ 1 x Nses

            arguments
                this mlvg.Inspector
                data {mustBeNumeric}
                pet_fqfn {mustBeText}
            end

            Nses = length(pet_fqfn);
            dstruct.wb = nan(1, Nses);
            dstruct.gm = nan(1, Nses);
            dstruct.wm = nan(1, Nses);
            dstruct.brainstem = nan(1, Nses);
            dstruct.cerebellum = nan(1, Nses);
            dstruct.subcortex = nan(1, Nses);

            % weight by anatomy
            for ifqfn = 1:length(pet_fqfn)
                try
                    w = this.weights_all_schaefer(pet_fqfn(ifqfn));  % ~ 309 x 1

                    w_wb = w .* ascol(ismember(this.indices_schaef, this.indices_wb));
                    w_wb = w_wb / sum(w_wb, "all");
                    w_wb = this.symmetrize_cortex(w_wb);
                    dstruct.wb(ifqfn) = median(data(:, ifqfn), "omitnan", Weights=ascol(w_wb));

                    w_gm = w .* ascol(ismember(this.indices_schaef, this.indices_gm));
                    w_gm = w_gm / sum(w_gm, "all");
                    w_gm = this.symmetrize_cortex(w_gm);
                    dstruct.gm(ifqfn) = median(data(:, ifqfn), "omitnan", Weights=ascol(w_gm));

                    w_wm = w .* ascol(ismember(this.indices_schaef, this.indices_wm));
                    w_wm = w_wm / sum(w_wm, "all");
                    w_wm = this.symmetrize_cortex(w_wm);
                    dstruct.wm(ifqfn) = median(data(:, ifqfn), "omitnan", Weights=ascol(w_wm));

                    w_bs = w .* ascol(ismember(this.indices_schaef, this.indices_brainstem));
                    w_bs = w_bs / sum(w_bs, "all");
                    w_bs = this.symmetrize_cortex(w_bs);
                    dstruct.brainstem(ifqfn) = median(data(:, ifqfn), "omitnan", Weights=ascol(w_bs));

                    w_cbm = w .* ascol(ismember(this.indices_schaef, this.indices_cerebellum));
                    w_cbm = w_cbm / sum(w_cbm, "all");
                    w_cbm = this.symmetrize_cortex(w_cbm);
                    dstruct.cerebellum(ifqfn) = median(data(:, ifqfn), "omitnan", Weights=ascol(w_cbm));

                    w_sc = w .* ascol(ismember(this.indices_schaef, this.indices_subcortex));
                    w_sc = w_sc / sum(w_sc, "all");
                    w_sc = this.symmetrize_cortex(w_sc);
                    dstruct.subcortex(ifqfn) = median(data(:, ifqfn), "omitnan", Weights=ascol(w_sc));
                catch ME
                    handwarning(ME)
                end
            end            
        end
    end

    methods (Static)

        function inspect_martinv1(martinv1, opts)
            arguments
                martinv1 {mustBeText}  % array of fqfn
                opts.plot_style = "raincloud"
            end

            this = mlvg.Inspector();
            data = this.load(martinv1);

            %% heatmap

            if strcmpi(opts.plot_style, "heatmap")
                % Sort sensors by median or mean to highlight the clustering
                medians = median(data, 2);
                [~, sortIdx] = sort(medians, 'descend');
                data_sorted = data(sortIdx, :);

                % Create heatmap
                figure('Position', [100, 100, 800, 1200]);
                imagesc(data_sorted');
                colorbar;
                xlabel('Schaefer parcel (sorted)');
                ylabel('v_1 (mL/cm^3)');
                title('v_1 Heatmap');
            end

            %% raincloud

            if strcmp(opts.plot_style, "raincloud")
                M = size(data, 1);
                figure('Position', [100, 100, 1200, 1400]);

                % Plot quantiles as thin lines
                subplot(1, 2, 1);
                quantiles = prctile(data, [10 25 50 75 90], 2);
                hold on;
                for m = 1:M
                    color = this.parc_color(m, M);

                    % Plot median as thicker line
                    plot(quantiles(m, 3), m, 'o', 'Color', color, 'MarkerSize', 2);
                    % Plot IQR
                    plot([quantiles(m, 2), quantiles(m, 4)], [m, m], '-', ...
                        'Color', color, 'LineWidth', 1.5);
                    % Plot 10-90 percentile
                    plot([quantiles(m, 1), quantiles(m, 5)], [m, m], '-', ...
                        'Color', color, 'LineWidth', 0.5);
                end
                ylim([0.5, M+0.5]);  % Add 0.5 padding to match imagesc behavior
                set(gca, 'YDir', 'reverse');
                ylabel('Schaefer parcel');
                xlabel('v_1 (mL/cm^3)');
                title('v_1 quantiles');

                % Add a compact heatmap showing density
                subplot(1, 2, 2);
                imagesc(data);
                colormap(gca, viridis);  % Apply viridis only to current axes
                ylabel('Schaefer parcel');
                xlabel('[^{15}O]CO PET session')
                c = colorbar;
                clim([dipmin(data), min(dipmax(data), 0.5)]);
                ylabel(c, 'v_1 (mL/cm^3)', 'Rotation', 270);
                title('Raw Data Heatmap');

                fontsize(scale=1.6)
            end
        end

        function inspect_huang1980(fdg, opts)
            arguments
                fdg {mustBeText}  % array of fqfn
                opts.plot_style = "raincloud"
                opts.measure {mustBeText} = "K_i"
            end

            this = mlvg.Inspector();

            tracer_name = "[^{18}F]FDG";
            switch convertStringsToChars(opts.measure)
                case "k_1"
                    measure_index = [];
                    measure_name = "k_1";
                    measure_units = "1/min";
                    converter = @(x) x * 60;
                case "k_1/k_2"
                    measure_index = 1;
                    measure_name = "k_1/k_2";
                    measure_units = "min/min";
                    converter = @(x) x;
                case "k_2"
                    measure_index = 2;
                    measure_name = "k_2";
                    measure_units = "1/min";
                    converter = @(x) x * 60;
                case "k_3"
                    measure_index = [];
                    measure_name = "k_3";
                    measure_units = "1/min";
                    converter = @(x) x * 60;
                case "k_3/k_4"
                    measure_index = 3;
                    measure_name = "k_3/k_4";
                    measure_units = "min/min";
                    converter = @(x) x;
                case "k_4"
                    measure_index = 4;
                    measure_name = "k_4";
                    measure_units = "1/min";
                    converter = @(x) x * 60;
                case "t_0"
                    measure_index = 5;
                    measure_name = "t_0";
                    measure_units = "s";
                    converter = @(x) x;
                case "sigma"
                    measure_index = 6;
                    measure_name = "\sigma";
                    measure_units = "Bq/Bq";
                    converter = @(x) x;
                case "K_i"
                    measure_index = [];
                    measure_name = "K_i";
                    measure_units = "mL hg^{-1}min^{-1}";  % v_1 ~ mL/cm^3
                    converter = @(x) x * 6000 / this.BRAIN_DENSITY;
                case "CMR_{glc}"
                    measure_index = [];
                    measure_name = "CMR_{glc}";
                    measure_units = "\mu mol hg^{-1}min^{-1}";
                    converter = @(x) x * 6000 / this.BRAIN_DENSITY;
                otherwise
                    error("mlvg:ValueError", stackstr());
            end

            if strcmpi(opts.measure, "k_1")
                [data_k1_k2,found] = this.load(fdg, measure_index=1);
                data_k2 = this.load(fdg, measure_index=2, converter=converter);
                data = data_k1_k2 .* data_k2;
            elseif strcmpi(opts.measure, "k_3")
                [data_k3_k4,found] = this.load(fdg, measure_index=3);
                data_k4 = this.load(fdg, measure_index=4, converter=converter);
                data = data_k3_k4 .* data_k4;
            elseif strcmpi(opts.measure, "K_i")
                [data_k1_k2,found] = this.load(fdg, measure_index=1);
                data_k2 = this.load(fdg, measure_index=2);
                data_k1 = data_k1_k2 .* data_k2;
                data_k3_k4 = this.load(fdg, measure_index=3);
                data_k4 = this.load(fdg, measure_index=4);
                data_k3 = data_k3_k4 .* data_k4;
                %data_k3 = data_k3_k4 .* 0.0063;  % Table 1.  Huang 1980
                data_v1 = this.load_matching_v1(fdg);
                assert(length(data_v1) == length(data_k1))

                denom = data_k2 + data_k3;
                data = data_v1 .* data_k1 .* data_k3 ./ denom;
                data = converter(data);
            elseif strcmpi(opts.measure, "CMR_{glc}")
                [data_k1_k2,found] = this.load(fdg, measure_index=1);
                data_k2 = this.load(fdg, measure_index=2);
                data_k1 = data_k1_k2 .* data_k2;
                data_k3_k4 = this.load(fdg, measure_index=3);
                data_k4 = this.load(fdg, measure_index=4);
                data_k3 = data_k3_k4 .* data_k4;
                data_v1 = this.load_matching_v1(fdg);
                assert(length(data_v1) == length(data_k1))

                denom = data_k2 + data_k3;
                data_Ki = data_v1 .* data_k1 .* data_k3 ./ denom;
                glc = this.load_matching_glc(fdg);
                data = data_Ki .* asrow(glc) / this.LC;
                data = converter(data);
            else
                [data,found] = this.load(fdg, measure_index=measure_index, converter=converter);
            end

            subs = string();
            fdg_found = fdg(found);
            for fidx = 1:length(fdg_found)
                f = fdg_found(fidx);
                s = mlpipeline.Bids.filename2struct(f);
                subs(fidx) = s.sub;
            end
            assert(size(subs, 2) == size(data, 2))

            %% raincloud

            if strcmp(opts.plot_style, "raincloud")
                M = size(data, 1);
                figure('Position', [100, 100, 2265, 1400]);

                % Plot quantiles as thin lines
                subplot(1, 3, 1);
                quantiles = prctile(data, [10 25 50 75 90], 2);
                hold on;
                for m = 1:M
                    color = this.parc_color(m, M);

                    % Plot median as thicker line
                    plot(quantiles(m, 3), m, 'o', 'Color', color, 'MarkerSize', 2);
                    % Plot IQR
                    plot([quantiles(m, 2), quantiles(m, 4)], [m, m], '-', ...
                        'Color', color, 'LineWidth', 1.5);
                    % Plot 10-90 percentile
                    plot([quantiles(m, 1), quantiles(m, 5)], [m, m], '-', ...
                        'Color', color, 'LineWidth', 0.5);
                end
                % if strcmpi(opts.measure, "k_3")
                %     xlim([0, 20 / 60])
                % end
                % if strcmpi(opts.measure, "k_3/k_4")
                %     xlim([0, 50 / 60])
                % end
                ylim([0.5, M+0.5]);  % Add 0.5 padding to match imagesc behavior
                set(gca, 'YDir', 'reverse');
                ylabel('Schaefer parcel');
                xlabel(sprintf('%s (%s)', measure_name, measure_units));
                title(sprintf('%s quantiles', measure_name));

                % Add a compact heatmap showing density
                h = subplot(1, 3, [2,3]);
                imagesc(data);                
                xticks(1:length(subs));  % Set the x-axis ticks to correspond to each column                
                xticklabels(subs);  % Label each column with strings from your array
                colormap(h, viridis);  % Apply viridis only to current axes
                ylabel('Schaefer parcel');
                xlabel(sprintf('%s PET session', tracer_name))
                c = colorbar;
                if strcmpi(opts.measure, "K_i")
                    clim([dipmin(data), 30]);
                else
                    clim([dipmin(data), dipmax(data)]);
                end
                ylabel(c, sprintf('%s (%s)', measure_name, measure_units), 'Rotation', 270);
                title('Raw Data Heatmap');

                fontsize(scale=1.6)
                set(h, 'XTickLabelRotation', 45, 'FontSize', 8);
            end
        end

        function inspect_mintun1984(oo, opts)
            arguments
                oo {mustBeText}  % array of fqfn
                opts.plot_style = "raincloud"
                opts.measure {mustBeText} = "oef"
                opts.fix_v_post logical = false
            end

            this = mlvg.Inspector();

            tracer_name = "[^{15}O]O_2";
            switch convertStringsToChars(opts.measure)
                case "oef"
                    measure_index = 1;
                    measure_name = "OEF";
                    measure_units = "Bq/Bq";
                    converter = @(x) x;
                case "f_wom"
                    measure_index = 2;
                    measure_name = "f_{wom} at 90 s";
                    measure_units = "mL/mL";
                    converter = @(x) x;
                case "v_post"
                    measure_index = 3;
                    measure_name = "(v_{post} + 0.5 v_{cap}) / v_1";
                    measure_units = "mL/mL";
                    converter = @(x) x;
                case "t_0"
                    measure_index = 4;
                    measure_name = "t_0";
                    measure_units = "s";
                    converter = @(x) x;
                case "tau_d"
                    measure_index = 5;
                    measure_name = "\tau_d";
                    measure_units = "s";
                    converter = @(x) x;
                case "sigma"
                    measure_index = 6;
                    measure_name = "\sigma";
                    measure_units = "Bq/Bq";
                    converter = @(x) x;
                case "CMRO_2"
                    measure_index = [];
                    measure_name = "CMRO_2";
                    measure_units = "\mu mol hg^{-1}min^{-1}";
                    converter = @(x) x * 6000 / this.BRAIN_DENSITY;
                otherwise
                    error("mlvg:ValueError", stackstr());
            end

            if strcmpi(opts.measure, "oef") && opts.fix_v_post
                data_oef = this.load(oo, measure_index=1, converter=converter);
                data_v_post = this.load(oo, measure_index=3, converter=converter);
                data = data_oef .* data_v_post / 0.835;
            elseif strcmpi(opts.measure, "CMRO_2")
                data_oef = this.load(oo, measure_index=1);
                data_v_post = this.load(oo, measure_index=3);
                data_f = this.load_matching_f(oo);
                assert(length(data_f) == length(data_oef))

                o2_content = this.load_matching_o2_content(oo);
                data = data_f .* data_oef .* (data_v_post / 0.835) .* o2_content;
                data = converter(data);
            else
                data = this.load(oo, measure_index=measure_index, converter=converter);
            end

            %% raincloud

            if strcmp(opts.plot_style, "raincloud")
                M = size(data, 1);
                figure('Position', [100, 100, 1200, 1400]);

                % Plot quantiles as thin lines
                subplot(1, 2, 1);
                quantiles = prctile(data, [10 25 50 75 90], 2);
                hold on;
                for m = 1:M
                    color = this.parc_color(m, M);

                    % Plot median as thicker line
                    plot(quantiles(m, 3), m, 'o', 'Color', color, 'MarkerSize', 2);
                    % Plot IQR
                    plot([quantiles(m, 2), quantiles(m, 4)], [m, m], '-', ...
                        'Color', color, 'LineWidth', 1.5);
                    % Plot 10-90 percentile
                    plot([quantiles(m, 1), quantiles(m, 5)], [m, m], '-', ...
                        'Color', color, 'LineWidth', 0.5);
                end
                ylim([0.5, M+0.5]);  % Add 0.5 padding to match imagesc behavior
                set(gca, 'YDir', 'reverse');
                ylabel('Schaefer parcel');
                xlabel(sprintf('%s (%s)', measure_name, measure_units));
                title(sprintf('%s quantiles', measure_name));

                % Add a compact heatmap showing density
                subplot(1, 2, 2);
                imagesc(data);
                colormap(gca, viridis);  % Apply viridis only to current axes
                ylabel('Schaefer parcel');
                xlabel(sprintf('%s PET session', tracer_name))
                c = colorbar;
                clim([dipmin(data), dipmax(data)]);
                ylabel(c, sprintf('%s (%s)', measure_name, measure_units), 'Rotation', 270);
                title('Raw Data Heatmap');

                fontsize(scale=1.6)
            end
        end

        function inspect_raichle1983(ho, opts)
            arguments
                ho {mustBeText}  % array of fqfn
                opts.plot_style = "raincloud"
                opts.measure {mustBeText} = "f"
            end

            tracer_name = "[^{15}O]H_2O";
            switch lower(convertStringsToChars(opts.measure))
                case "f"
                    measure_index = 1;
                    measure_name = "f";
                    measure_units = "mL/min/cm^3";
                    converter = @(x) 60*x;
                case "lambda"
                    measure_index = 2;
                    measure_name = "\lambda";
                    measure_units = "cm^3/mL";
                    converter = @(x) x;
                case "ps"
                    measure_index = 3;
                    measure_name = "PS";
                    measure_units = "mL/min/cm^3";
                    converter = @(x) 60*x;
                case "t_0"
                    measure_index = 4;
                    measure_name = "t_0";
                    measure_units = "s";
                    converter = @(x) x;
                case "sigma"
                    measure_index = 5;
                    measure_name = "\sigma";
                    measure_units = "Bq/Bq";
                    converter = @(x) x;
                otherwise
                    error("mlvg:ValueError", stackstr());
            end

            this = mlvg.Inspector();
            data = this.load(ho, measure_index=measure_index, converter=converter);

            %% heatmap

            if strcmpi(opts.plot_style, "heatmap")
                % Sort sensors by median or mean to highlight the clustering
                medians = median(data, 2);
                [~, sortIdx] = sort(medians, 'descend');
                data_sorted = data(sortIdx, :);

                % Create heatmap
                figure('Position', [100, 100, 800, 1200]);
                imagesc(data_sorted');
                colorbar;
                xlabel('Schaefer parcel (sorted)');
                ylabel(sprintf('%s (%s)', measure_name, measure_units));
                title(sprintf('%s Heatmap', measure_name));
            end

            %% raincloud

            if strcmp(opts.plot_style, "raincloud")
                M = size(data, 1);
                figure('Position', [100, 100, 1200, 1400]);

                % Plot quantiles as thin lines
                subplot(1, 2, 1);
                quantiles = prctile(data, [10 25 50 75 90], 2);
                hold on;
                for m = 1:M
                    color = this.parc_color(m, M);

                    % Plot median as thicker line
                    plot(quantiles(m, 3), m, 'o', 'Color', color, 'MarkerSize', 2);
                    % Plot IQR
                    plot([quantiles(m, 2), quantiles(m, 4)], [m, m], '-', ...
                        'Color', color, 'LineWidth', 1.5);
                    % Plot 10-90 percentile
                    plot([quantiles(m, 1), quantiles(m, 5)], [m, m], '-', ...
                        'Color', color, 'LineWidth', 0.5);
                end
                ylim([0.5, M+0.5]);  % Add 0.5 padding to match imagesc behavior
                set(gca, 'YDir', 'reverse');
                ylabel('Schaefer parcel');
                xlabel(sprintf('%s (%s)', measure_name, measure_units));
                title(sprintf('%s quantiles', measure_name));

                % Add a compact heatmap showing density
                subplot(1, 2, 2);
                imagesc(data);
                colormap(gca, viridis);  % Apply viridis only to current axes
                ylabel('Schaefer parcel');
                xlabel(sprintf('%s PET session', tracer_name))
                c = colorbar;
                clim([dipmin(data), dipmax(data)]);
                ylabel(c, sprintf('%s (%s)', measure_name, measure_units), 'Rotation', 270);
                title('Raw Data Heatmap');

                fontsize(scale=1.6)
            end
        end
    end

    %% PRIVATE

    properties
        indices_schaef_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
