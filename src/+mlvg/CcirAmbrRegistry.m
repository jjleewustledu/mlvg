classdef (Sealed) CcirAmbrRegistry < handle & mlpipeline.StudyRegistry
    %% CCIRAMBRREGISTRY
    %  
    %  Created 05-Feb-2023 00:55:13 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.13.0.2126072 (R2022b) Update 3 for MACI64.  Copyright 2023 John J. Lee.
        
    methods (Static)
        function this = instance()
            persistent uniqueInstance      
            if (isempty(uniqueInstance))
                this = mlvg.CcirAmbrRegistry();
                uniqueInstance = this;
            else
                this = uniqueInstance;
            end
        end
    end

    %% PRIVATE

    methods (Access = private)
        function this = CcirAmbrRegistry()
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
