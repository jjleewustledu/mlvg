classdef QuadraticAerobicGlycolysisKit < handle & mlpet.AbstractAerobicGlycolysisKit
    %% QUADRATICAEROBICGLYCOLYSISKIT is a factory implementing quadratic parameterization of kinetic rates using
    %  emissions.  See also papers by Videen, Herscovitch.  This implementation supports CCIR_01211.
    %  
    %  Created 25-Jan-2022 20:01:06 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.11.0.1837725 (R2021b) Update 2 for MACI64.  Copyright 2022 John J. Lee.
    
	methods (Static)
        function construct(varargin)
            %% CONSTRUCT
            %  e.g.:  construct('cbv', 'subjectsExpr', 'sub-S58163*', 'Nthreads', 1, 'region', 'wholebrain', 'aifMethods', 'idif')
            %  e.g.:  construct('cbv', 'debug', true)
            %  @param required physiolog is char, e.g., cbv, cbf, cmro2, cmrglc.
            %  @param subjectsExpr is text, e.g., 'sub-S58163*'.
            %  @param region is char, e.g., wholebrain, voxels.
            %  @param debug is logical.
            %  @param Nthreads is numeric|char.
            
            import mlvg.*
            import mlvg.QuadraticAerobicGlycolysisKit.*

            % global
            registry = MatlabRegistry.instance(); %#ok<NASGU>
            setenv('SUBJECTS_DIR', fullfile(getenv('SINGULARITY_HOME'), 'CCIR_012211', 'derivatives'))
            setenv('PROJECTS_DIR', getenv('SINGULARITY_HOME'))
            setenv('DEBUG', '')
            setenv('NOPLOT', '')
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'physiology', @istext)
            addParameter(ip, 'subjectsExpr', 'sub-108293', @istext)
            addParameter(ip, 'region', 'voxels', @istext)
            addParameter(ip, 'debug', ~isempty(getenv('DEBUG')), @islogical)
            addParameter(ip, 'Nthreads', 1, @(x) isnumeric(x) || istext(x))
            parse(ip, varargin{:})
            ipr = ip.Results;
            if istext(ipr.Nthreads)
                ipr.Nthreads = str2double(ipr.Nthreads);
            end 
            
            % switch strategy
            switch ipr.physiology
                case 'cbv'
                    tracer = 'oc';
                    metric = 'vs';
                    region = ipr.region;
                    construction = @QuadraticAerobicGlycolysisKit.constructCbvByRegion;                            
                case 'cbf'
                    tracer = 'ho';
                    metric = 'fs';
                    region = ipr.region;
                    construction = @QuadraticAerobicGlycolysisKit.constructCbfByRegion;
                case 'cmro2'
                    tracer = 'oo';
                    metric = 'os';
                    region = ipr.region;
                    construction = @QuadraticAerobicGlycolysisKit.constructCmro2ByRegion;
                otherwise
                    error('mlvg:RuntimeError', 'QuadraticAerobicGlycolysisKit.construct.ipr.physiology->%s', ipr.physiology)
            end
            
            % construct            
            pwd1 = pushd(getenv('SUBJECTS_DIR'));
            mlvg.Ccir1211Registry.instance('initialize')
            theSessionData = QuadraticAerobicGlycolysisKit.constructSessionData( ...
                metric, ...
                'subjectsExpr', ipr.subjectsExpr, ...
                'tracer', tracer, ...
                'debug', ipr.debug, ...
                'region', region); % length(theSessionData) ~ 60
            if ipr.Nthreads > 1                
                parfor (p = 1:length(theSessionData), ipr.Nthreads)
                    try
                        construction(theSessionData(p)); %#ok<PFBNS>
                    catch ME
                        handwarning(ME)
                    end
                end
            elseif ipr.Nthreads == 1
                for p = length(theSessionData):-1:1
                    try
                        construction(theSessionData(p)); % RAM ~ 3.3 GB
                    catch ME
                        handwarning(ME)
                    end
                end
            end            
            popd(pwd1);
        end
        function constructCbfByRegion(varargin)
            %% CONSTRUCTCBFBYREGION
            %  @param required sessionData is mlpipeline.ISessionData.
            %  @return cbf on filesystem.
            
            this = mlvg.QuadraticAerobicGlycolysisKit(varargin{:});            
            Region = [upper(this.sessionData.region(1)) this.sessionData.region(2:end)];

            pwd0 = pushd(this.sessionData.subjectPath);
            fs_ = this.(['buildFsBy' Region])();             
            cbf_ = this.fs2cbf(fs_);
            cbf_.save() % save ImagingContext2            
            popd(pwd0);
        end
        function constructCbvByRegion(varargin)
            %% CONSTRUCTCBVBYREGION
            %  @param required sessionData is mlpipeline.ISessionData.
            %  @return cbv on filesystem.
            
            this = mlvg.QuadraticAerobicGlycolysisKit(varargin{:});
            Region = [upper(this.sessionData.region(1)) this.sessionData.region(2:end)];

            pwd0 = pushd(this.sessionData.subjectPath);            
            vs_ = this.(['buildVsBy' Region])(); 
            cbv_ = this.vs2cbv(vs_);
            cbv_.save() % save ImagingContext2
            popd(pwd0);
        end 
        function constructCmro2ByRegion(varargin)
            %% CONSTRUCTCMRO2BYREGION
            %  @param required sessionData is mlpipeline.ISessionData.
            %  @return cmro2 on filesystem.
            %  @return oef on filesystem.
            
            this = mlvg.QuadraticAerobicGlycolysisKit(varargin{:});            
            this.constructPhysiologyDateOnly('cbf', ...
                'subjectFolder', this.sessionData.subjectFolder, ...
                'region', this.sessionData.region, ...
                'sessionData', this.sessionData)
            this.constructPhysiologyDateOnly('cbv', ...
                'subjectFolder', this.sessionData.subjectFolder, ...
                'region', this.sessionData.region, ...
                'sessionData', this.sessionData)
            Region = [upper(this.sessionData.region(1)) this.sessionData.region(2:end)];

            pwd0 = pushd(this.sessionData.subjectPath);             
            os_ = this.(['buildOsBy' Region])();   
            cbf_ = this.sessionData.cbfOnAtlas( ...
                'typ', 'mlfourd.ImagingContext2', ...
                'dateonly', true, ...
                'tags', [this.blurTag this.sessionData.regionTag]);
            [cmro2_,oef_] = this.os2cmro2(os_, cbf_, this.model);
            cmro2_.save() % save ImagingContext2
            oef_.save()            
            popd(pwd0);
        end  
        function ic = constructPhysiologyDateOnly(varargin)
            %% e.g., constructPhysiologyDateOnly('cbv', 'sessionData', obj) % pwd includes sub-108293            
            %  e.g., constructPhysiologyDateOnly('cbv', 'sessionData', obj, ...
            %                                    'workpath', '/path/to/sub-108293/resampling_restricted')
            %  e.g., constructPhysiologyDateOnly('cbv', 'sessionData', obj, ...
            %                                    'workpath', '/path/to/sub-108293/resampling_restricted', ...
            %                                    'filepatt', 'ocdt2021*_111_voxels.4dfp.hdr')
            
            import mlvg.QuadraticAerobicGlycolysisKit
            
            reg = mlvg.Ccir1211Registry.instance();
            ip = inputParser;
            addRequired(ip, 'physiology', @istext)
            addParameter(ip, 'sessionData', [], @(x) isa(x, 'mlpipeline.ISessionData'))
            addParameter(ip, 'workpath', pwd, @isfolder)
            addParameter(ip, 'subjectFolder', '', @istext) 
            addParameter(ip, 'filepatt', '', @istext)
            addParameter(ip, 'atlTag', reg.atlasTag, @istext)
            addParameter(ip, 'blurTag', reg.blurTag, @istext)
            addParameter(ip, 'region', 'voxels', @istext)
            parse(ip, varargin{:})
            ipr = ip.Results;
            if isempty(ipr.subjectFolder)
                ipr.subjectFolder = reg.workpath2subject(ipr.workpath);
            end
            if isempty(ipr.filepatt)

                ipr.filepatt = sprintf('%sdt*%s%s_%s.4dfp.hdr', ipr.physiology, ipr.atlTag, ipr.blurTag, ipr.region);
                %ipr.filepatt = sprintf('%s_ses-*_trac-*_proc-%s_pet%s_%s.4dfp.hdr', ...
                %    ipr.subjectFolder, ipr.physiology, ipr.atlTag, ipr.region);
            end
            
            pwd0 = pushd(ipr.workpath);
            g = globT(ipr.filepatt);
            if isempty(g); return; end
            
            %% segregate by dates
            
            m = containers.Map;            
            for ig = 1:length(g)
                if contains(g{ig}, ipr.sessionData.defects)
                    continue
                end
                dstr = AbstractAerobicGlycolysisKit.physiologyObjToDatetimeStr(g{ig}, 'dateonly', true);
                if ~lstrfind(m.keys, dstr)
                    m(dstr) = g(ig); % cell
                else
                    m(dstr) = [m(dstr) g{ig}];
                end
            end 
            
            %% average scans by dates
            
            for k = asrow(m.keys)
                fns = m(k{1});
                ic = mlfourd.ImagingContext2(fns{1});
                ic = ic.zeros();
                icfp = strrep(ic.fileprefix, ...
                    QuadraticAerobicGlycolysisKit.physiologyObjToDatetimeStr(fns{1}), ...
                    QuadraticAerobicGlycolysisKit.physiologyObjToDatetimeStr(fns{1}, 'dateonly', true));
                if isfile([icfp '.4dfp.hdr'])
                    continue
                end
                ic_count = 0;
                for fn = fns
                    incr = mlfourd.ImagingContext2(fn{1});
                    if dipsum(incr) > 0
                        ic = ic + incr;
                        ic_count = ic_count + 1;
                    end
                end
                ic = ic / ic_count;
                ic.fileprefix = icfp;
                ic.save()
            end
            
            %%
            
            popd(pwd0);
        end  
        function theSD = constructSessionData(varargin)
            
            import mlvg.*
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'metric', @istext)
            addParameter(ip, 'subjectsExpr', 'sub-*', @istext)
            addParameter(ip, 'sessionsExpr', 'ses-*', @istext)
            addParameter(ip, 'tracer', '', @istext)
            addParameter(ip, 'debug', ~isempty(getenv('DEBUG')), @islogical)
            addParameter(ip, 'region', 'voxels', @istext)
            addParameter(ip, 'scanIndex', 1:2, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            idx = 1;
            subsPath = getenv('SUBJECTS_DIR');
            pwd1 = pushd(subsPath);
            subjects = globFoldersT(ipr.subjectsExpr); % e.g., 'sub-*'
            for sub = subjects
                pwd0 = pushd(fullfile(subsPath, sub{1}, 'pet'));
                subd = SubjectData('subjectFolder', sub{1});
                sesfs = globFoldersT(ipr.sessionsExpr); % e.g., 'ses-*'

                for ses = sesfs
                    for scan_idx = ipr.scanIndex
                        try
                            sesd = SessionData( ...
                                'studyData', StudyData(), ...
                                'projectData', ProjectData('sessionStr', ses{1}), ...
                                'subjectData', subd, ...
                                'sessionFolder', ses{1}, ...
                                'scanIndex', scan_idx, ...
                                'tracer', upper(ipr.tracer), ...
                                'ac', true, ...
                                'region', ipr.region, ...
                                'metric', ipr.metric);            
                            if ~isfile(sesd.wmparc1OnAtlas)
                                mlpet.AbstractAerobicGlycolysisKit.constructWmparc1OnAtlas(sesd);
                            end
                            tracerfn = sesd.([lower(sesd.tracer) 'OnAtlas']);
                            if ~isfile(tracerfn)
                                sesd.jitOnAtlas(tracerfn)
                            end
                            theSD(idx) = sesd; %#ok<AGROW>
                            idx = idx + 1;
                        catch ME
                            if strcmpi('mlvg:ValueError:getScanFolder', ME.identifier)
                                continue
                            end
                            handwarning(ME)
                        end
                    end
                end
                popd(pwd0);
            end
            popd(pwd1);
        end
        function dt = physiologyObjToDatetime(obj)
            ic = mlfourd.ImagingContext2(obj);            
            ss = split(ic.fileprefix, '_');
            re = regexp(ss{2}, 'ses-(?<datetime>\d{14})\w*', 'names');
            dt = datetime(re.datetime, 'InputFormat', 'yyyyMMddHHmmss');
        end
        function dtstr = physiologyObjToDatetimeStr(varargin)
            import mlvg.QuadraticAerobicGlycolysisKit 
            ip = inputParser;
            addRequired(ip, 'obj', @(x) ~isempty(x))
            addParameter(ip, 'dateonly', false, @islogical)
            parse(ip, varargin{:})
            ipr = ip.Results;  
            if ipr.dateonly
                dtstr = [datestr(QuadraticAerobicGlycolysisKit.physiologyObjToDatetime(ipr.obj), 'yyyymmdd') '000000'];
            else
                dtstr = datestr(QuadraticAerobicGlycolysisKit.physiologyObjToDatetime(ipr.obj), 'yyyymmddHHMMSS') ;
            end
        end
    end

    properties 
        aifMethods
        indexCliff
        model
        sessionData
    end
    
    properties (Dependent)
        atlasTag
        blurTag
        dataFolder % e.g., resampling_restricted
        dataPath
        registry
        regionTag
        subjectFolder % e.g., sub-108293
        subjectPath
        tags
    end 

    methods
        
        %% GET
        
        function g = get.atlasTag(~)
            g = mlvg.Ccir1211Registry.instance.atlasTag;
        end
        function g = get.blurTag(~)
            g = mlvg.Ccir1211Registry.instance.blurTag;
        end
        function g = get.dataFolder(this)
            g = this.sessionData.dataFolder;
        end  
        function g = get.dataPath(this)
            g = this.sessionData.dataPath;
        end  
        function g = get.registry(this)
            g = this.sessionData.registry;
        end
        function g = get.regionTag(this)
            g = this.sessionData.regionTag;
        end
        function g = get.subjectFolder(this)
            g = this.sessionData.subjectFolder;
        end
        function g = get.subjectPath(this)
            g = this.sessionData.subjectPath;
        end
        function g = get.tags(this)
            g = strcat(this.blurTag, this.regionTag);
        end
        
        %%

        function fs_ = buildFsByVoxels(this, varargin)
            %% BUILDFSBYVOXELS
            %  @return fs in R^4 as mlfourd.ImagingContext2, without saving to filesystems.  
            
            import mloxygen.QuadraticNumericRaichle1983
            
            ensuredir(this.dataPath);
            pwd0 = pushd(this.dataPath);  
                                    
            brain = this.sessionData.brainOnAtlas('typ', 'mlfourd.ImagingContext2'); 
            devkit = mlpet.ScannerKit.createFromSession(this.sessionData);             
            scanner = devkit.buildScannerDevice();
            scannerBrain = scanner.volumeAveraged(brain.binarized());
            arterial = this.buildAif(devkit, scanner, scannerBrain);
            
            fs_ = copy(brain.fourdfp);
            fs_.filepath = this.dataPath;
            fs_.fileprefix = this.fsOnAtlas('typ', 'fp', 'tags', this.tags);

            % solve Raichle
            fprintf('%s\n', datestr(now))
            fprintf('starting mlraichle.QuadraticAerobicGlycolysisKit.buildFsByVoxels\n')
            raichle = QuadraticNumericRaichle1983.createFromDeviceKit( ...
                devkit, ...
                'scanner', scanner, ...
                'arterial', arterial, ...
                'roi', brain.binarized());  
            raichle = raichle.solve();
            this.model = raichle;

            % insert Raichle solutions into fs
            fs_.img = raichle.fs('typ', 'single');
                
            fs_ = mlfourd.ImagingContext2(fs_);
            popd(pwd0);
        end 
        function os_ = buildOsByVoxels(this, varargin)
            %% BUILDOSBYVOXELS
            %  @return os in R^4 as mlfourd.ImagingContext2, without saving to filesystems.  
                    
            import mloxygen.QuadraticNumericMintun1984
            
            ensuredir(this.dataPath);
            pwd0 = pushd(this.dataPath);  
                                    
            brain = this.sessionData.brainOnAtlas('typ', 'mlfourd.ImagingContext2'); 
            devkit = mlpet.ScannerKit.createFromSession(this.sessionData);            
            scanner = devkit.buildScannerDevice(); 
            scannerBrain = scanner.volumeAveraged(brain.binarized()); 
            arterial = this.buildAif(devkit, scanner, scannerBrain);
            
            os_ = copy(brain.fourdfp);
            os_.filepath = this.dataPath;
            os_.fileprefix = this.osOnAtlas('typ', 'fp', 'tags', this.tags);

            % solve Mintun
            fprintf('%s\n', datestr(now))
            fprintf('starting mlraichle.QuadraticAerobicGlycolysisKit.buildOsByVoxels\n')
            mintun = QuadraticNumericMintun1984.createFromDeviceKit( ...
                devkit, ...
                'scanner', scanner, ...
                'arterial', arterial, ...
                'roi', brain.binarized());  
            mintun = mintun.solve();
            this.model = mintun;

            % insert Raichle solutions into fs
            os_.img = mintun.os('typ', 'single');

            os_ = mlfourd.ImagingContext2(os_);
            popd(pwd0);
        end
        function vs_ = buildVsByVoxels(this, varargin)
            %% BUILDVSBYVOXELS
            %  @return v1_ in R^ as mlfourd.ImagingContext2, without saving to filesystems.  
            %  @return cbv_ in R^3, without saving.
            
            import mloxygen.QuadraticNumericMartin1987
            
            ensuredir(this.dataPath);
            pwd0 = pushd(this.dataPath);                                    
            
            brain = this.sessionData.brainOnAtlas('typ', 'mlfourd.ImagingContext2');
            devkit = mlpet.ScannerKit.createFromSession(this.sessionData);             
            scanner = devkit.buildScannerDevice();
            scannerBrain = scanner.volumeAveraged(brain.binarized());
            arterial = this.buildAif(devkit, scanner, scannerBrain);
            
            vs_ = copy(brain.fourdfp);
            vs_.filepath = this.dataPath;
            vs_.fileprefix = this.vsOnAtlas('typ', 'fp', 'tags', this.tags);

            % solve Martin
            fprintf('%s\n', datestr(now))
            fprintf('starting mlraiche.QuadraticAerobicGlycolysisKit.buildVsByVoxels\n')
            martin = QuadraticNumericMartin1987.createFromDeviceKit( ...
                devkit, ...
                'scanner', scanner, ...
                'arterial', arterial, ...
                'roi', brain.binarized());
            martin = martin.solve();  
            this.model = martin;

            % insert Martin solutions into fs
            vs_.img = martin.vs('typ', 'single');

            vs_ = mlfourd.ImagingContext2(vs_);
            popd(pwd0);
        end        
        function obj = metricOnAtlas(this, metric, varargin)
            %% METRICONATLAS 
            %  @param required metric is char.
            %  @param datetime is datetime or char, .e.g., '20200101000000' | ''.
            %  @param dateonly is logical.
            %  @param tags is char, e.g., 'b43_wmparc1', default ''.
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'metric', @ischar)
            addParameter(ip, 'datetime', this.sessionData.datetime, @(x) isdatetime(x) || ischar(x))
            addParameter(ip, 'dateonly', false, @islogical)
            addParameter(ip, 'tags', '', @istext)
            parse(ip, metric, varargin{:})
            ipr = ip.Results;
            ipr.metric = lower(ipr.metric);

            if ~isempty(ipr.tags)
                ipr.tags = strcat("_", strip(ipr.tags, "_"));
            end
            if ischar(ipr.datetime)
                adatestr = ipr.datetime;
            end
            if isdatetime(ipr.datetime)
                if ipr.dateonly
                    adatestr = [datestr(ipr.datetime, 'yyyymmdd') '000000'];
                else
                    adatestr = datestr(ipr.datetime, 'yyyymmddHHMMSS');
                end
            end
            
            % e.g., derivatives/sub-108293/pet/sub-108293_ses-20210421144815_trc-oc_proc-dyn_pet_on_T1w.nii.gz
            fqfn = fullfile( ...
                this.dataPath, ...
                sprintf('%s_ses-%s_trc-%s_proc-%s_pet_%s%s%s', ...
                        this.subjectFolder, ...
                        adatestr, ...
                        this.metric2tracer(ipr.metric), ...
                        ipr.metric, ...
                        this.atlasTag, ...
                        ipr.tags, ...
                        this.sessionData.filetypeExt));
            obj  = this.sessionData.fqfilenameObject(fqfn, varargin{:});
        end	

        function this = QuadraticAerobicGlycolysisKit(varargin)
            %% QUADRATICAEROBICGLYCOLYSISKIT 
            %  Args:
            %      sessionData (ISessionData): session-specific objects
            %      aifMethods: containers.Map
            
            this = this@mlpet.AbstractAerobicGlycolysisKit(varargin{:})
            
            am = containers.Map;
            am('CO') = 'idif';
            am('OC') = 'idif';
            am('OO') = 'idif';
            am('HO') = 'idif';
            am('FDG') = 'idif';

            ip = inputParser;
            ip.KeepUnmatched = true;
            ip.PartialMatching = false;
            addRequired(ip, 'sessionData', @(x) isa(x, 'mlpipeline.ISessionData'))
            addParameter(ip, 'indexCliff', [], @isnumeric)
            addParameter(ip, 'aifMethods', am, @(x) isa(x, 'containers.Map'))
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.sessionData = ipr.sessionData;
            this.indexCliff = ipr.indexCliff;
            this.aifMethods = ipr.aifMethods;
            
            this.resetModelSampler()
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
