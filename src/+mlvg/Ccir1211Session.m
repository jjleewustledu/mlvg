classdef Ccir1211Session < mlpipeline.SessionData2 & handle
    %% line1
    %  line2
    %  
    %  Created 06-Feb-2023 21:31:45 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.13.0.2126072 (R2022b) Update 3 for MACI64.  Copyright 2023 John J. Lee.
        
    properties
        defects = {}
    end

    methods
        function this = Ccir1211Session(varargin)
            this = this@mlpipeline.SessionData2(varargin{:});
            this.registry_ = mlsiemens.VisionRegistry.instance();
        end
    end

    %% PROTECTED

    methods (Access = ?mlpipeline.SessionData2)
        function buildRadmeasurements(this)
            this.radMeasurements_ = mlpet.CCIRRadMeasurements2.createFromSession( ...
                this.mediator_);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
