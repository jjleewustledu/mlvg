classdef (Sealed) Ccir1211Registry < handle & mlpipeline.StudyRegistry
	%% CCIR1211REGISTRY 

	%  $Revision$
 	%  was created 15-Oct-2015 16:31:41
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlvg/src/+mlvg.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.
 	    
    properties
        atlasTag = 'on_T1w'
        blurTag = ''
        comments = ''
        Ddatetime0 % seconds
        dicomExtension = '.dcm'
        ignoredExperiments = {}
        noclobber = true
        normalizationFactor = 1
        numberNodes
        projectFolder = 'CCIR_01211';
        referenceTracer = 'FDG'
        scatterFraction = 0
        T = 0 % sec at the start of artery_interpolated used for model but not described by scanner frames
        tracerList = {'oc' 'oo' 'ho' 'fdg'}
        umapType = 'ct'
        stableToInterpolation = true
        voxelTime = 60 % sec
        wallClockLimit = 168*3600 % sec
    end
    
    properties (Dependent)
        earliestCalibrationDatetime        
        projectsDir
        rawdataDir
        subjectsDir
        subjectsJson
        tBuffer
    end
    
    methods % GET        
        function g = get.earliestCalibrationDatetime(~)
            %g = datetime(2015,1,1, 'TimeZone', 'local'); % accomodates sub-S33789
            g = datetime(2016,7,19, 'TimeZone', 'local');
        end
        function g = get.projectsDir(~)
            g = getenv('SINGULARITY_HOME');
        end 
        function x = get.rawdataDir(this)
            x = fullfile(this.projectsDir, this.projectFolder, 'rawdata');
        end
        function g = get.subjectsDir(this)
            g = fullfile(this.projectsDir, this.projectFolder, 'derivatives');
        end    
        function g = get.subjectsJson(this)
            if isempty(this.subjectsJson_)
                this.subjectsJson_ = jsondecode( ...
                    fileread(fullfile(this.projectsDir, this.projectFolder, 'constructed_20210225.json')));
            end
            g = this.subjectsJson_;
        end
        function g = get.tBuffer(this)
            g = max(0, -this.Ddatetime0) + this.T;
        end
    end

    methods

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
        function pth = sourcedataPathToDerivativesPath(pth)
            pth = strrep(pth, 'sourcedata', 'derivatives');
        end
        function sub = workpath2subject(pth)
            c = fileparts2cell(pth);
            msk = contains(c, 'sub-');
            sub = c(msk);
        end
    end  

    %% PRIVATE

    properties (Access = private)
        subjectsJson_
    end
    
	methods (Access = private)		  
 		function this = Ccir1211Registry()
 		end
    end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

