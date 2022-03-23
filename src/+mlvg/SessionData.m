classdef SessionData < mlnipet.MetabolicSessionData
	%% SESSIONDATA  

	%  $Revision$
 	%  was created 15-Feb-2016 01:51:37
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlvg/src/+mlvg.
 	%% It was developed on Matlab 9.0.0.307022 (R2016a) Prerelease for MACI64.    
        
    methods (Static)
        function t = consoleTaus(tracer)
            reg = mlvg.Ccir1211Registry.instance();
            t = reg.consoleTaus(tracer);
        end 
        function this = create(varargin)
            %  Args:
            %      folders (text): <proj folder>/derivatives/<sub folder>/<ses folder>/pet[/file], 
            %                      e.g., 'CCIR_01211/derivatives/sub-108293/ses-20210421/pet'
            %                      e.g., '$SINGULARITY_HOME/CCIR_01211/derivatives/sub-108293/ses-20210421/pet/sub-108293_ses-20210421155709_trc-fdg_proc-dyn_pet_on_T1w.nii.gz'
            %      ignoreFinishMark (logical): default := false
            %      reconstructionMethod (text): e.g., 'e7', 'NiftyPET'
            
            ip = inputParser;
            addRequired(ip, 'folders', @istext);
            addParameter(ip, 'ignoreFinishMark', false, @islogical);
            addParameter(ip, 'reconstructionMethod', 'e7', @istext);
            addParameter(ip, 'studyRegistry', mlvg.Ccir1211Registry.instance());
            parse(ip, varargin{:});
            [ipr,b,ic] = adjustIpr(ip.Results);
    
            this = mlvg.SessionData( ...
                'studyData', ipr.studyRegistry, ...
                'projectData', mlvg.ProjectData('projectFolder', ipr.prjfold), ...
                'subjectData', mlvg.SubjectData('subjectFolder', ipr.subfold), ...
                'sessionFolder', ipr.sesfold, ...
                'scanFolder', ipr.scnfold, ...
                'bids', b, ...
                'imagingContext', ic);
            this.ignoreFinishMark = ipr.ignoreFinishMark;  
            this.reconstructionMethod = ipr.reconstructionMethod;
            
            function [ipr,b,ic] = adjustIpr(ipr)
                ss = strsplit(ipr.folders, filesep);  
                ipr.prjfold = '';
                ipr.subfold = '';
                ipr.sesfold = '';
                ipr.scnfold = '';
                if any(contains(ss, 'CCIR_'))
                    ipr.prjfold = ss{contains(ss, 'CCIR')}; % 1st match
                end
                if any(contains(ss, 'sub-'))
                    ipr.subfold = ss{contains(ss, 'sub-')};
                end
                if any(contains(ss, 'ses-'))
                    ipr.sesfold = ss{contains(ss, 'ses-')};
                end
                if any(contains(ss, '-Converted-'))
                    ipr.scnfold = ss{contains(ss, '-Converted-')};
                end

                b = []; ic = [];
                if isfolder(ipr.folders)
                    b = mlvg.Ccir1211Bids('destinationPath', ipr.folders);
                end
                if isfile(ipr.folders)
                    b = mlvg.Ccir1211Bids('destinationPath', myfileparts(ipr.folders));
                    ic = mlfourd.ImagingContext2(ipr.folders);
                end
            end
        end
        function obj = struct2sessionData(obj) %#ok<INUSD> 
            error('mlvg:NotImplementedError', 'SessionData.struct2sessionData')
        end
    end

    properties (Constant)
        STUDY_CENSUS_XLSX_FN = ''
    end
    
    properties
        defects = {}
        tracers = {'fdg' 'ho' 'oo' 'oc'}
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

        dataPath
        dataFolder

        bids
        imagingContext
        registry
    end

    methods
             
        %% GET/SET

        function g    = get.projectsDir(this)
            if ~isempty(this.studyData_)
                g = this.studyData_.projectsDir;
                return
            end
            if ~isempty(this.projectData_)
                g = this.projectData_.projectsDir;
                return
            end
            error('mlvg:RuntimeError', 'SessionData.get.projectsDir')
        end
        function this = set.projectsDir(this, s)
            assert(istext(s))
            if ~isempty(this.studyData_)
                this.studyData_.projectsDir = s;
                return
            end
            if ~isempty(this.projectData_)
                this.projectData_.projectsDir = s;
                return
            end
            error('mlvg:RuntimeError', 'SessionData.get.projectsDir')
        end
        function g    = get.projectsPath(this)
            g = this.projectsDir;
        end
        function g    = get.projectsFolder(this)
            g = mybasename(this.projectsDir);
        end     
        function g    = get.projectPath(this)
            g = fullfile(this.projectsPath, this.projectFolder);
        end
        function g    = get.projectFolder(this)
            g = this.projectData.projectFolder;
        end
        
        function g    = get.subjectsDir(this)
            g = fullfile(this.projectPath, 'derivatives', '');
        end
        function this = set.subjectsDir(this, s)
            assert(istext(s));
            this.projectData_.projectsDir = myfileparts(s);
        end
        function g    = get.subjectsPath(this)
            g = this.subjectsDir;
        end
        function g    = get.subjectsFolder(this)
            g = mybasename(this.subjectsDir);
        end 
        function g    = get.subjectPath(this)
            g = fullfile(this.subjectsDir, this.subjectFolder);
        end
        function g    = get.subjectFolder(this)
            g = this.subjectData.subjectFolder;
        end        

        function g    = get.sessionsDir(this)
            g = this.subjectPath;
        end
        function g    = get.sessionsPath(this)
            g = this.sessionsDir;
        end
        function g    = get.sessionsFolder(this)
            g = mybasename(this.sessionsDir);
        end
        function g    = get.sessionPath(this)
            g = fullfile(this.sessionsDir, this.sessionFolder);
        end
        function this = set.sessionPath(this, s)
            assert(istext(s));
            [pth,this.sessionFolder] = myfileparts(s);
            [pth_,this.subjectData.subjectFolder] = myfileparts(pth);
            pth__ = myfileparts(pth_); % drop derivatives
            [this.projectsDir,this.projectData.projectFolder] = myfileparts(pth__);
        end
        function g    = get.sessionFolder(this)
            g = this.sessionFolder_;
        end        
        function this = set.sessionFolder(this, s)
            assert(istext(s));
            this.sessionFolder_ = s;            
        end    
        
        function g    = get.scansDir(this)
            g = this.sessionPath;
        end
        function g    = get.scansPath(this)
            g = this.scansDir;
        end
        function g    = get.scansFolder(this)
            g = this.sessionFolder;
        end
        function g    = get.scanPath(this)
            g = fullfile(this.scansPath, this.scanFolder);
        end
        function this = set.scanPath(this, s)
            assert(istext(s));
            [this.sessionPath,this.scanFolder] = myfileparts(s);
        end
        function g    = get.scanFolder(this)
            if (~isempty(this.scanFolder_))
                g = this.scanFolder_;
                return
            end

            %% KLUDGE for bootstrapping
            
            if isempty(this.tracer_) || isempty(this.attenuationCorrected_)
                g = '';
                dt = datetime(datestr(now));
                for globbed = globFoldersT(fullfile(this.sessionPath, '*_DT*.000000-Converted-*'))
                    base = mybasename(globbed{1});
                    re = regexp(base, ...
                        '\S+_DT(?<yyyy>\d{4})(?<mm>\d{2})(?<dd>\d{2})(?<HH>\d{2})(?<MM>\d{2})(?<SS>\d{2})\.\d{6}-Converted\S*', ...
                        'names');
                    assert(~isempty(re))
                    dt1 = datetime(str2double(re.yyyy), str2double(re.mm), str2double(re.dd), ...
                        str2double(re.HH), str2double(re.MM), str2double(re.SS));
                    if dt1 < dt
                        dt = dt1; % find earliest scan
                        g = base;
                    end                    
                end                
                return
            end
            dtt = mlpet.DirToolTracer( ...
                'tracer', fullfile(this.sessionPath, this.tracer_), ...
                'ac', this.attenuationCorrected_);            
            assert(~isempty(dtt.dns));
            try
                g = dtt.dns{this.scanIndex};
            catch ME
                if length(dtt.dns) < this.scanIndex 
                    error('mlnipet:ValueError:getScanFolder', ...
                        'SessionData.getScanFolder().this.scanIndex->%s', mat2str(this.scanIndex))
                else
                    rethrow(ME)
                end
            end
        end
        function this = set.scanFolder(this, s)
            assert(istext(s))
            this = this.setScanFolder(s);
        end

        function g    = get.dataPath(this)
            g = fullfile(this.subjectPath, this.dataFolder, '');
        end
        function g    = get.dataFolder(~)
            g = 'resampling_restricted';
        end

        function g    = get.bids(this)
            g = copy(this.bids_);
        end
        function g    = get.imagingContext(this)
            g = copy(this.imagingContext_);
        end
        function g    = get.registry(this)
            g = this.registry_;
        end

        %%
        
        function t = alternativeTaus(this)
            t = mlvg.SessionData.consoleTaus(this.tracer);
        end
        function p = petPointSpread(~, varargin)
            vis = mlsiemens.VisionRegistry.instance();
            p = vis.petPointSpread(varargin{:});
        end
        function tracerRawdataLocation(~)
            error('mlvg:NotImplementedError', 'SessionData.tracerRawdataLocation');
        end
        
      	function this = SessionData(varargin)
 			this = this@mlnipet.MetabolicSessionData(varargin{:}); 
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'bids', []);
            addParameter(ip, 'imagingContext', []);
            addParameter(ip, 'registry', mlvg.Ccir1211Registry.instance());
            parse(ip, varargin{:});   
            ipr = ip.Results;
            this.bids_ = ipr.bids;
            this.imagingContext_ = ipr.imagingContext;
            this.registry_ = ipr.registry;
            if isempty(this.tracer_) && ~isempty(this.bids_) && ~isempty(this.imagingContext_)
                this.tracer_ = this.bids_.obj2tracer(this.imagingContext_);
            end

            this.ReferenceTracer = 'FDG';
            if isempty(this.studyData_)
                this.studyData_ = mlvg.StudyData();
            end
            if isempty(this.projectData_)
                this.projectData_ = mlvg.ProjectData('sessionStr', this.sessionFolder);
            end
        end
    end

    %% PRIVATE

    properties (Access = private)
        bids_
        imagingContext_
        registry_
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

