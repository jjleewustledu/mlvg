classdef (Sealed) Ccir1211Registry < handle & mlpipeline.IStudyRegistry
	%% STUDYREGISTRY 

	%  $Revision$
 	%  was created 15-Oct-2015 16:31:41
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlvg/src/+mlvg.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.
 	    
    methods (Static)
        function t = consoleTaus(tracer)
            %% see also t0_and_dt()
                      
            switch (upper(tracer))
                case 'FDG'
                    t = [5*ones(1,24) 20*ones(1,9) 60*ones(1,10) 300*ones(1,9)];
                case {'OO' 'HO'}
                    t = [3*ones(1,23) 5*ones(1,6) 10*ones(1,8) 30*ones(1,6)];
                case {'OC' 'CO'}
                    t = [15 60*ones(1,5)];
                otherwise
                    error('mlvg:IndexError', 'SessionData.consoleTaus.tracer->%s', tracer);
            end
        end 
        function this = instance(varargin)
            %% INSTANCE
            %  @param optional qualifier is char \in {'initialize' ''}
            
            ip = inputParser;
            addOptional(ip, 'qualifier', '', @ischar)
            parse(ip, varargin{:})
            
            persistent uniqueInstance
            if (strcmp(ip.Results.qualifier, 'initialize'))
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

    properties (Constant)
        PREFERRED_TIMEZONE = 'America/Chicago'
    end

    properties
        atlasTag = '111'
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
    
    methods
        
        %% GET
        
        function g = get.earliestCalibrationDatetime(~)
            %g = datetime(2015,1,1, 'TimeZone', 'America/Chicago'); % accomodates sub-S33789
            g = datetime(2016,7,19, 'TimeZone', 'America/Chicago');
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

        %%  

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

