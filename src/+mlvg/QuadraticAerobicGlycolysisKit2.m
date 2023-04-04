classdef QuadraticAerobicGlycolysisKit2 < handle & mlpet.AbstractAerobicGlycolysisKit2
    %% QUADRATICAEROBICGLYCOLYSISKIT2 is a factory implementing quadratic parameterization of kinetic rates using
    %  emissions.  See also papers by Videen, Herscovitch.  This implementation supports CCIR_01211.
    %  
    %  Created 04-Apr-2023 12:25:22 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.14.0.2206163 (R2023a) for MACI64.  Copyright 2023 John J. Lee.
   

    %% PROTECTED

    methods (Access = {?mlpet.AbstractAerobicGlycolysisKit2, ?mlvg.QuadraticAerobicGlycolysisKit2})
        function this = QuadraticAerobicGlycolysisKit2(varargin)
            %% QUADRATICAEROBICGLYCOLYSISKIT 
            %  Args:
            %      immediator (ImagingMediator): session-specific objects
            %      aifMethods: containers.Map
            
            this = this@mlpet.AbstractAerobicGlycolysisKit2(varargin{:})
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
