classdef Ccir1211Mediator < handle & mlpipeline.ImagingMediator
    %% CCIR1211MEDIATOR provides a mediator design pattern for project CCIR1211 (cf. GoG pp. 276-278).  
    %  As a mediator, it separates and manages data-conceptual entities previously squashed into class 
    %  hierarchies such as that for mlvg.SessionData.
    %
    %  It also provides a prototype design pattern for use by abstract factories like mlkinetics.BidsKit 
    %  (cf. GoF pp. 90-91, 117).  For prototypes, call initialize(obj) using obj understood by 
    %  mlfourd.ImagingContext2.  Delegates data-conceptual functionality to mlvg.{Ccir1211Scan, Ccir1211Session, 
    %  Ccir1211Subject, Ccir1211Project, Ccir1211Study}.
    %  
    %  Created 06-Feb-2023 19:33:37 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.13.0.2126072 (R2022b) Update 3 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        defects
        metric
    end

    methods
        function this = Ccir1211Mediator(varargin)
            %% Args must be understandable by mlfourd.ImagingContext2.

            this = this@mlpipeline.ImagingMediator(varargin{:});
            try
                this.initialize();
            catch ME
                handwarning(ME)
            end
        end
        function [this,T] = findProximal(this, index, patt)
            arguments
                this mlvg.Ccir1211Mediator
                index double = 1
                patt {mustBeTextScalar} = "sub-*_ses-*_trc-fdg_proc-dyn_pet.nii.gz"
            end

            % table:  dt, fqfn
            g = glob(fullfile(this.subjectsPath, "sub-*", "ses-*", this.scanFolder, patt));
            T = cell2table(cellfun(@(x) this.datetime_bids_filename(x), g, UniformOutput=false));
            T = addvars(T, g);
            T.Properties.VariableNames = ["dt", "fqfn"];

            % table:  dt, sep, fqfn; sorted by sep
            current = contains(T.fqfn, this.imagingContext.filename);
            sep = days(abs(T.dt - T.dt(current)));
            T = addvars(T, sep, After="dt", NewVariableNames="sep");
            T = sortrows(T, "sep"); % ascending

            % find proximal Ccir1211Mediator
            fqfn1 = T.fqfn{1+index};
            this = mlvg.Ccir1211Mediator(fqfn1);
        end
        function this = initialize(this, varargin)
            this.buildImaging(varargin{:});
            
            this.bids_ = mlvg.Ccir1211Bids( ...
                originationPath=this.scanPath, ...
                projectPath=this.projectPath, ...
                subjectFolder=this.subjectFolder);
            this.imagingContext_ = this.ensureFiniteImagingContext(this.imagingContext_);
            this.imagingAtlas_ = this.bids_.atlas_ic;
            try
                this.imagingDlicv_ = this.bids_.dlicv_ic;
            catch ME
                handwarning(ME)
            end          
        end
        function icd = prepare_derivatives(this, ic)
            arguments
                this mlvg.Ccir1211Mediator
                ic {mustBeNonempty} = this.imagingContext_
            end

            icd = this.bids.prepare_derivatives(ic);
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
                imcontext = this.imagingContext_
            end
            if ~isempty(imcontext)
                this.imagingContext_ = imcontext;
            end

            this.scanData_ = mlvg.Ccir1211Scan(this, dataPath=this.imagingContext_.filepath);
            this.sessionData_ = mlvg.Ccir1211Session(this, dataPath=this.scansPath);
            this.subjectData_ = mlvg.Ccir1211Subject(this, dataPath=this.sessionsPath);
            this.projectData_ = mlvg.Ccir1211Project(this, dataPath= ...            
                this.omit_bids_folders(this.subjectsPath));
            this.studyData_ = mlvg.Ccir1211Study(this, mlvg.Ccir1211Registry.instance());

            % additional assembly
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
