classdef Ccir1211Scan < mlpipeline.ScanData2 & handle
    %% line1
    %  line2
    %  
    %  Created 06-Feb-2023 21:32:20 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.13.0.2126072 (R2022b) Update 3 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        defects = {}
    end

    properties (Dependent)
        isotope
        reconstructionMethod
        tracer % lower-case
    end

    methods % GET
        function g = get.isotope(this)
            switch lower(this.tracer)
                case 'fdg'
                    g = '18F';
                case {'co' 'ho' 'oc' 'oh' 'oo'}
                    g = '15O';
                otherwise
                    error("mlvg:ValueError", stackstr())
            end
        end
        function g = get.tracer(this)
            ic = this.mediator_.imagingContext;
            if ~contains(ic.fileprefix, "pet")
                g = "none";
                return
            end
            re = regexp(ic.fileprefix, "_trc-(?<trc>\w+)_", "names");
            g = lower(re.trc);
        end
        function g = get.reconstructionMethod(this)
            g = 'e7';
        end
    end

    methods
        function this = Ccir1211Scan(varargin)
            this = this@mlpipeline.ScanData2(varargin{:});
        end
        function dt = datetime(this, varargin)
            dt = this.datetime_bids_filename(varargin{:});
            deltadt = seconds(this.mediator_.timeOffsetConsole);
            dt = dt - deltadt;
        end
        function ic = metricOnAtlas(this, metric, tags)
            %% METRICONATLAS forms an ImagingContext2 with modality->metric
            %  and adding tags and atlasTag.

            arguments
                this mlvg.Ccir1211Scan
                metric {mustBeTextScalar} = 'unknown'
                tags {mustBeTextScalar} = ''
            end

            s = this.mediator_.bids.filename2struct(this.mediator_.imagingContext.fqfn);
            s.modal = metric;
            s.tag = strcat(s.tag, tags);
            fqfn = fullfile(this.mediator_.scanPath, this.mediator_.bids.struct2filename(s));
            ic  = mlfourd.ImagingContext2(fqfn);
        end
        function t = taus(this, trc)
            arguments
                this mlvg.Ccir1211Scan
                trc {mustBeTextScalar} = this.tracer
            end
            t = this.consoleTaus(trc);
        end
    end

    methods (Static)
        function t = consoleTaus(tracer)
            %% entered into console protocol
            %  see also t0_and_dt()
                      
            switch lower(tracer)
                case 'fdg'
                    t = [5*ones(1,24) 20*ones(1,9) 60*ones(1,10) 300*ones(1,9)];
                case {'oo' 'ho'}
                    t = [3*ones(1,23) 5*ones(1,6) 10*ones(1,8) 30*ones(1,6)];
                case {'oc' 'co'}
                    t = [15 60*ones(1,5)];
                otherwise
                    error('mlvg:IndexError', 'SessionData.consoleTaus.tracer->%s', tracer);
            end
        end 
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
