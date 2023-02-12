classdef Ccir1211Subject < mlpipeline.SubjectData2 & handle
    %% line1
    %  line2
    %  
    %  Created 06-Feb-2023 21:31:12 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.13.0.2126072 (R2022b) Update 3 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        defects = {}
    end

    methods
        function this = Ccir1211Subject(varargin)
            this = this@mlpipeline.SubjectData2(varargin{:});
            setenv('SUBJECTS_DIR', this.subjectsDir);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
