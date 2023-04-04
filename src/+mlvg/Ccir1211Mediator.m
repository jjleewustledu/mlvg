classdef Ccir1211Mediator < handle & mlpipeline.ImagingMediator
    %% CCIR1211MEDIATOR provides a mediator design pattern for project CCIR1211.
    %  
    %  Created 06-Feb-2023 19:33:37 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.13.0.2126072 (R2022b) Update 3 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        defects
        metric
    end

    properties (Constant)
        BIDS_FOLDERS = {'derivatives', 'rawdata', 'sourcedata'};
    end

    properties (Dependent)
        anatPath
        derivAnatPath
        derivativesPath
        derivPetPath
        mriPath
        petPath
        sourcedataPath
        sourceAnatPath
        sourcePetPath

        atlas_ic
        dlicv_ic
        flair_ic
        T1_ic % FreeSurfer
        T1_on_t1w_ic
        t1w_ic
        t2w_ic
        tof_ic
        tof_on_t1w_ic
        wmparc_ic % FreeSurfer
        wmparc_on_t1w_ic % FreeSurfer
    end

    methods % GET/SET
        function g = get.anatPath(this)
            g = this.bids.anatPath;
        end      
        function g = get.derivAnatPath(this)
            g = this.bids.derivAnatPath;
        end
        function g = get.derivativesPath(this)
            g = this.bids.derivativesPath;
        end
        function g = get.derivPetPath(this)
            g = this.bids.derivPetPath;
        end        
        function g = get.mriPath(this)
            g = this.bids.mriPath;
        end
        function g = get.petPath(this)
            g = this.bids.petPath;
        end
        function g = get.sourcedataPath(this)
            g = this.bids.sourcedataPath;
        end
        function g = get.sourceAnatPath(this)
            g = this.bids.sourceAnatPath;
        end
        function g = get.sourcePetPath(this)
            g = this.bids.sourcePetPath;
        end
        
        function g = get.atlas_ic(this)
            g = this.bids.atlas_ic;
        end  
        function g = get.dlicv_ic(this)
            g = this.bids.dlicv_ic;
        end
        function g = get.flair_ic(this)
            g = this.bids.flair_ic;
        end
        function g = get.T1_ic(this) % FreeSurfer
            g = this.bids.T1_ic;
        end
        function g = get.T1_on_t1w_ic(this)
            g = this.bids.T1__on_t1w_ic;
        end
        function g = get.t1w_ic(this)
            g = this.bids.t1w_ic;
        end
        function g = get.t2w_ic(this)
            g = this.bids.t2w_ic;
        end
        function g = get.tof_ic(this)
            g = this.bids.tof_ic;
        end
        function g = get.tof_on_t1w_ic(this)
            g = this.bids.tof_on_t1w_ic;
        end
        function g = get.wmparc_ic(this) % FreeSurfer
            g = this.bids.wmparc_ic;
        end
        function g = get.wmparc_on_t1w_ic(this)
            g = this.bids.wmparc_on_t1w_ic;
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
        function icd = prepare_derivatives(this, ic)
            icd = this.bids.prepare_derivatives(ic);
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
