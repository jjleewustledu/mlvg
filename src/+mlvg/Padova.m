classdef Padova < handle
    %% line1
    %  line2
    %  
    %  Created 09-Mar-2022 00:22:11 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.11.0.1873467 (R2021b) Update 3 for MACI64.  Copyright 2022 John J. Lee.
    
    properties
        dcvmap
        ralmap
    end

    methods
        function m = bloodSuckerMap(~)
            m = containers.Map;
            m('k1') = struct('min',  1,     'max',  4,     'init',  3,    'sigma', 0.05); % alpha
            m('k2') = struct('min',  0.05,  'max',  0.5,   'init',  0.2,  'sigma', 0.05); % beta in 1/sec
            m('k3') = struct('min',  1,     'max',   3,    'init',  1.5,  'sigma', 0.05); % p
            m('k4') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % dp2 for 2nd bolus
            m('k5') = struct('min',  0,     'max', 100,    'init', 10,    'sigma', 0.05); % t0 in sec
            m('k6') = struct('min',  0.95,  'max',   1,    'init',  0.995,'sigma', 0.05); % steady-state fraction in (0, 1), for rising baseline
            m('k7') = struct('min',  0.02,  'max',   0.15, 'init',  0.05, 'sigma', 0.05); % recirc fraction < 0.5, for 2nd bolus
            m('k8') = struct('min', 15,     'max',  90,    'init', 30,    'sigma', 0.05); % recirc delay in sec
            m('k9') = struct('min',  0,     'max',   0.1,  'init',  0.01, 'sigma', 0.05); % baseline amplitude fraction \approx 0.05
        end
        function histogram_and_plot(this)
            a = [];
            b = [];
            p = [];
            ssf = [];
            loss = [];
            for k = asrow(this.dcvmap.keys)
                ral = this.solve_ecat(this.dcvmap(k{1}));
                plot(ral)
                this.ralmap(k{1}) = ral;
                a = [a ral.strategy.ks(1)];
                b = [b ral.strategy.ks(2)];
                p = [p ral.strategy.ks(3)];
                ssf = [ssf ral.strategy.ks(6)];
                loss = [loss ral.strategy.loss];
            end
            figure; histogram(a); title('histogram alpha');
            figure; histogram(b); title('histogram beta');
            figure; histogram(p); title('histogram p');
            figure; histogram(ssf); title('histogram ssf');
            figure; histogram(loss); title('histogram loss');
            saveFigures(pwd, 'closeFigure', false)
        end
        function plot_ecats(this)
            for k = asrow(this.dcvmap.keys)
                ral = this.solve_ecat(this.dcvmap(k{1}));
                plot(ral)
                title(k, 'FontSize', 12)
            end
        end
        function ral = solve_ecat(this, wc)
            ral = mlswisstrace.RadialArteryLee2021( ...
                'tracer', 'HO', ...
                'kernel', 1, ...
                'model_kind', '2bolus', ...
                'map', this.bloodSuckerMap(), ...
                'Measurement', wc);
            ral = ral.solve();
            plot(ral)
        end

        function this = Padova(varargin)
            %% PADOVA 
            %  Args:
            %      arg1 (its_class): Description of arg1.
            
            ip = inputParser;
            addParameter(ip, "arg1", "", @istext)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            pad = load('/Users/jjlee/Library/CloudStorage/Box-Box/Tommaso - General/IDIF/OneDrive_1_3-7-2022/IDIF_H2O_all.mat');
            this.dcvmap = containers.Map('KeyType', 'double', 'ValueType', 'any');
            for ip = 1:length(pad.IDIF_all)
                t = pad.IDIF_all(ip).t_IDIF_emi;
                activity = pad.IDIF_all(ip).IDIF .* 2.^(-t/122.2416);
                activity = interp1(t, activity, t(1):t(end), 'linear', 0);
                activity = activity(1:60);
                this.dcvmap(ip) = asrow(activity);
            end
            this.ralmap = containers.Map('KeyType', 'double', 'ValueType', 'any');
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
