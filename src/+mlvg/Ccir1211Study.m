classdef Ccir1211Study < handle & mlpipeline.StudyData2
    %% line1
    %  line2
    %  
    %  Created 06-Feb-2023 21:28:27 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.13.0.2126072 (R2022b) Update 3 for MACI64.  Copyright 2023 John J. Lee.
    
    methods
        function this = Ccir1211Study(varargin)
            this = this@mlpipeline.StudyData2(varargin{:});
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
