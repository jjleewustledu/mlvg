classdef Ccir1211Mediator < handle & mlpipeline.ImagingMediator
    %% line1
    %  line2
    %  
    %  Created 06-Feb-2023 19:33:37 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.13.0.2126072 (R2022b) Update 3 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        defects
        metric
    end

    properties (Dependent)
        projectsDir % homolog of __Freesurfer__ subjectsDir
        projectsPath
        projectsFolder
        projectPath
        projectFolder % \in projectsFolder        
        
        subjectsDir % __Freesurfer__ convention
        subjectsPath 
        subjectsFolder 
        subjectPath
        subjectFolder % \in subjectsFolder  
        
        sessionsDir % __Freesurfer__ convention
        sessionsPath 
        sessionsFolder 
        sessionPath
        sessionFolder % \in projectFolder        
        
        scansDir % __Freesurfer__ convention
        scansPath 
        scansFolder 
        scanPath
        scanFolder % \in sessionFolder

        isotope
        reconstructionMethod
        tracer
    end

    methods % GET/SET
        function g = get.projectsDir(this)
            g = this.projectData_.projectsDir;
        end
        function g = get.projectsPath(this)
            g = this.projectData_.projectsPath;
        end
        function g = get.projectsFolder(this)
            g = this.projectData_.projectsFolder;
        end
        function g = get.projectPath(this)
            g = this.projectData_.projectPath;
        end
        function     set.projectPath(this, s)
            this.projectData_.projectPath = s;
        end
        function g = get.projectFolder(this)
            g = this.projectData_.projectFolder;
        end        
        function     set.projectFolder(this, s)
            this.projectData_.projectFolder = s;
        end   

        function g = get.subjectsDir(this)
            g = this.subjectData_.subjectsDir;
        end
        function g = get.subjectsPath(this)
            g = this.subjectData_.subjectsPath;
        end
        function g = get.subjectsFolder(this)
            g = this.subjectData_.subjectsFolder;
        end
        function g = get.subjectPath(this)
            g = this.subjectData_.subjectPath;
        end
        function     set.subjectPath(this, s)
            this.subjectData_.subjectPath = s;
        end
        function g = get.subjectFolder(this)
            g = this.subjectData_.subjectFolder;
        end        
        function     set.subjectFolder(this, s)
            this.subjectData_.subjectFolder = s;
        end   

        function g = get.sessionsDir(this)
            g = this.sessionData_.sessionsDir;
        end
        function g = get.sessionsPath(this)
            g = this.sessionData_.sessionsPath;
        end
        function g = get.sessionsFolder(this)
            g = this.sessionData_.sessionsFolder;
        end
        function g = get.sessionPath(this)
            g = this.sessionData_.sessionPath;
        end
        function     set.sessionPath(this, s)
            this.sessionData_.sessionPath = s;
        end
        function g = get.sessionFolder(this)
            g = this.sessionData_.sessionFolder;
        end        
        function     set.sessionFolder(this, s)
            this.sessionData_.sessionFolder = s;
        end

        function g = get.scansDir(this)
            g = this.scanData_.scansDir;
        end
        function g = get.scansPath(this)
            g = this.scanData_.scansPath;
        end
        function g = get.scansFolder(this)
            g = this.scanData_.scansFolder;
        end
        function g = get.scanPath(this)
            g = this.scanData_.scanPath;
        end
        function     set.scanPath(this, s)
            this.scanData_.scanPath = s;
        end
        function g = get.scanFolder(this)
            g = this.scanData_.scanFolder;
        end        
        function     set.scanFolder(this, s)
            this.scanData_.scanFolder = s;
        end

        function g = get.isotope(this)
            g = this.scanData_.isotope;
        end
        function g = get.reconstructionMethod(this)
            g = this.scanData_.reconstructionMethod;
        end
        function g = get.tracer(this)
            g = this.scanData_.tracer;
        end
    end

    methods
        function this = Ccir1211Mediator(varargin)
            %% Args must be understandable by mlfourd.ImagingContext2.

            this = this@mlpipeline.ImagingMediator(varargin{:});
            this.buildImaging();
            this.bids_ = mlvg.Ccir1211Bids( ...
                destinationPath=this.scanPath, ...
                projectPath=this.projectPath, ...
                subjectFolder=this.subjectFolder);
            this.imagingAtlas_ = this.bids_.atlas_ic;
            this.imagingDlicv_ = this.bids_.dlicv_ic;
        end
        function dt = datetime(this, varargin)
            dt = this.scanData_.datetime(varargin{:});
        end
        function imagingChanged(this, imdata)
            %% subclasses override to affect mlpipeline.ImagingData
            %  complexity of mediator design patterns arise here

            arguments
                this mlpipeline.ImagingMediator
                imdata mlpipeline.ImagingData
            end

            if imdata == this.scanData_
            elseif imdata == this.sessionData_
            elseif imdata == this.subjectData_
            elseif imdata == this.projectData_
            elseif imdata == this.studyData_
            else
                error('mlpipeline:ValueError', stackstr());
            end
        end
        function ic = metricOnAtlas(this, varargin)
            ic = this.scanData_.metricOnAtlas(varargin{:});
        end
        function ps = petPointSpread(this, varargin)
            reg = this.sessionData_.registry;
            ps = reg.petPointSpread(varargin{:});
        end
        function t = taus(this, varargin)
            t = this.scanData_.taus(varargin{:});
        end
        function ic = tracerOnAtlas(this, varargin)
            if endsWith(this.imagingContext.fileprefix, this.registry.atlasTag)
                ic = this.imagingContext;
                return
            end
            s = this.bids.filename2struct(this.imagingContext.fqfn);
            s.tag = this.atlasTag;
            fqfn = this.bids.struct2filename(s);
            ic = mlfourd.ImagingContext2(fqfn);
        end
    end

    methods (Static)
        function this = create(varargin)
            this = mlvg.Ccir1211Mediator(varargin{:});
        end
    end

    %% PROTECTED

    methods (Access = protected)        
        function buildImaging(this, imcontext)
            arguments
                this mlvg.Ccir1211Mediator
                imcontext = []
            end
            if ~isempty(imcontext)
                this.imagingContext_ = mlfourd.ImagingContext2(imcontext);
            end

            this.scanData_ = mlvg.Ccir1211Scan(this, dataPath=this.imagingContext_.filepath);
            this.sessionData_ = mlvg.Ccir1211Session(this, dataPath=this.scansPath);
            this.subjectData_ = mlvg.Ccir1211Subject(this, dataPath=this.sessionsPath);
            this.projectData_ = mlvg.Ccir1211Project(this, dataPath= ...            
                this.omit_bids_folders(this.subjectsPath));
            this.studyData_ = mlvg.Ccir1211Study(this, mlvg.Ccir1211Registry.instance());

            % additional assembly required?

        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
