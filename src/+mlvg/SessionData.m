classdef SessionData < mlnipet.MetabolicSessionData
	%% SESSIONDATA  

	%  $Revision$
 	%  was created 15-Feb-2016 01:51:37
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlvg/src/+mlvg.
 	%% It was developed on Matlab 9.0.0.307022 (R2016a) Prerelease for MACI64.    
    
    properties (Constant)
        STUDY_CENSUS_XLSX_FN = ''
    end
    
    properties
        registry
        tracers = {'fdg' 'ho' 'oo' 'oc'}
    end
    
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
        function this = create(varargin)
            % @param folders ~ <project folder>/<session folder>/<scan folder>, in getenv('SINGULARITY_HOME')
            % @param ignoreFinishMark is logical, default := false
            
            ip = inputParser;
            addRequired(ip, 'folders', @(x) isfolder(fullfile(getenv('SINGULARITY_HOME'), x)))
            addParameter(ip, 'ignoreFinishMark', false, @islogical);
            addParameter(ip, 'reconstructionMethod', 'e7', @ischar);
            parse(ip, varargin{:});
            ipr = adjustIpr(ip.Results);
    
            this = mlvg.SessionData( ...
                'studyData', mlvg.Ccir1211Registry.instance(), ...
                'projectData', mlvg.ProjectData('projectFolder', ipr.prjfold), ...
                'subjectData', mlvg.SubjectData('subjectFolder', ipr.subfold), ...
                'sessionFolder', ipr.sesfold, ...
                'scanFolder', ipr.scnfold);
            this.ignoreFinishMark = ipr.ignoreFinishMark;  
            this.reconstructionMethod = ipr.reconstructionMethod;          
            
            function ipr = adjustIpr(ipr)
                ss = strsplit(ipr.folders, filesep);
                if lstrfind(ss{1}, 'subjects')
                    p = mlvg.ProjectData();
                    ipr.subfold = ss{2};
                    ipr.prjfold = p.session2project(ss{3});
                    ipr.sesfold = ss{3};
                    if length(ss) >= 3
                        ipr.scnfold = ss{4};
                    else
                        ipr.scnfold = '';
                    end
                    return
                end                
                
                ipr.prjfold = ss{1};
                ipr.subfold = mlvg.SubjectData().sesFolder2subFolder(ss{2});
                ipr.sesfold = ss{2};
                if length(ss) >= 3
                    ipr.scnfold = ss{3};
                else
                    ipr.scnfold = '';
                end
            end
        end
        function sessd = struct2sessionData(sessObj)
            if (isa(sessObj, 'mlvg.SessionData'))
                sessd = sessObj;
                return
            end
            
            import mlvg.*;
            assert(isfield(sessObj, 'projectFolder'));
            assert(isfield(sessObj, 'sessionFolder'));
            assert(isfield(sessObj, 'sessionDate'));
            studyd = StudyData;
            sessp = fullfile(studyd.projectsDir, sessObj.projectFolder, sessObj.sessionFolder, '');
            sessd = SessionData('studyData', studyd, ...
                                'sessionPath', sessp, ...
                                'tracer', 'FDG', ...
                                'ac', true, ...
                                'sessionDate', sessObj.sessionDate);  
            if ( isfield(sessObj, 'parcellation') && ...
                ~isempty(sessObj.parcellation))
                sessd.parcellation = sessObj.parcellation;
            end
            if ( isfield(sessObj, 'region') && ...
                ~isempty(sessObj.region))
                sessd.region = sessObj.region;
            end
        end
    end

    methods
                
        %%
        
        function getStudyCensus(this, ~)
            error('mlvg:NotImplementedError', 'SessionData.studyCensus');
        end            
        function tracerRawdataLocation(this, ~)
            error('mlvg:NotImplementedError', 'SessionData.tracerRawdataLocation');
        end
        
      	function this = SessionData(varargin)
 			this = this@mlnipet.MetabolicSessionData(varargin{:}); 
            if isempty(this.studyData_)
                this.studyData_ = mlvg.StudyData();
            end
            this.ReferenceTracer = 'FDG';
            if isempty(this.projectData_)
                this.projectData_ = mlvg.ProjectData('sessionStr', this.sessionFolder);
            end
            
            %% registry
            
            this.registry = mlvg.Ccir1211Registry.instance();
            
            %% taus
            
            if (~isempty(this.scanFolder_) && isfile(this.jsonFilename, 'file'))
                j = jsondecode(fileread(this.jsonFilename));
                this.taus_ = j.taus';
            end
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

