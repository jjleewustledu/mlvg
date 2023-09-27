classdef Ccir1211Scan < handle & mlpipeline.ScanData2
    %% line1
    %  line2
    %  
    %  Created 06-Feb-2023 21:32:20 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.13.0.2126072 (R2022b) Update 3 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
    end

    methods
        function this = Ccir1211Scan(varargin)
            this = this@mlpipeline.ScanData2(varargin{:});
        end
    end

    methods (Static)
        function t = consoleTaus(tracer)
            %% entered into console protocol
            %  see also t0_and_dt()
                      
            if ~istext(tracer)
                t = 0;
                return
            end
            switch lower(tracer)
                case 'fdg'
                    t = [5*ones(1,24) 20*ones(1,9) 60*ones(1,10) 300*ones(1,9)];
                case {'oo' 'ho'}
                    t = [3*ones(1,23) 5*ones(1,6) 10*ones(1,8) 30*ones(1,6)];
                case {'oc' 'co'}
                    t = [15 60*ones(1,5)];
                otherwise
                    t = 0; % for non-PET imaging such as T1w
            end
        end 
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
 