classdef Ccir1211Json < mlpipeline.CcirJson
    %% CCIR1211JSON manages *.json and *.mat data repositories for CNDA-related data such as subjects, experiments, aliases.
    %  
    %  Created 23-Feb-2022 17:53:35 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.11.0.1873467 (R2021b) Update 3 for MACI64.  Copyright 2022 John J. Lee.
    
    methods (Static)        
        function S = loadConstructed()
            import mlvg.Ccir1211Json;
            S = jsondecode(fileread(fullfile(Ccir1211Json.projectPath, Ccir1211Json.filenameConstructed)));
        end
    end

    properties (Constant)
        filenameConstructed = 'constructed_20210225.json'      
        projectPath = fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211', '')
    end

    methods
        function this = Ccir1211Json(varargin)
            this = this@mlpipeline.CcirJson(varargin{:});
            %this.S_ = mlvg.Ccir1211Json.loadConstructed();   
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
