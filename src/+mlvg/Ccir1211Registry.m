classdef (Sealed) Ccir1211Registry < handle & mlpipeline.StudyRegistry
	%% CCIR1211REGISTRY 
    %
	%  $Revision$
 	%  was created 15-Oct-2015 16:31:41
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlvg/src/+mlvg.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.
 	    
    properties
        snakes
    end

    properties (Dependent)
        subjectsJson
    end
    
    methods % GET
        function g = get.subjectsJson(this)
            if isempty(this.subjectsJson_)
                this.subjectsJson_ = jsondecode( ...
                    fileread(fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211', 'constructed_20210225.json')));
            end
            g = this.subjectsJson_;
        end
    end
    
    methods (Static)
        function t = consoleTaus(tracer)
            t = mlvg.Ccir1211Scan.consoleTaus(tracer);
        end 
        function this = instance(reset)
            arguments
                reset = []
            end
            persistent uniqueInstance
            if ~isempty(reset)
                uniqueInstance = [];
            end
            if (isempty(uniqueInstance))
                this = mlvg.Ccir1211Registry();
                uniqueInstance = this;
            else
                this = uniqueInstance;
            end
        end        
    end  

    %% PRIVATE
    
	methods (Access = private)		  
 		function this = Ccir1211Registry()
            this.atlasTag = 'on_T1w';
            this.reconstructionMethod = 'e7';
            this.referenceTracer = 'FDG';
            this.tracerList = {'oc' 'oo' 'ho' 'fdg'};
            this.T = 0;
            this.umapType = 'ct';

            this.snakes.contractBias = 0.2;
            this.snakes.iterations = 80;
            this.snakes.smoothFactor = 0;
 		end
    end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

