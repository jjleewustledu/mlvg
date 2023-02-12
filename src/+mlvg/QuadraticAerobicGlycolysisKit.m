classdef QuadraticAerobicGlycolysisKit < handle & mlpet.AbstractAerobicGlycolysisKit2
    %% QUADRATICAEROBICGLYCOLYSISKIT is a factory implementing quadratic parameterization of kinetic rates using
    %  emissions.  See also papers by Videen, Herscovitch.  This implementation supports CCIR_01211.
    %  
    %  Created 25-Jan-2022 20:01:06 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.11.0.1837725 (R2021b) Update 2 for MACI64.  Copyright 2022 John J. Lee.
    
	methods (Static)
        function these = construct(varargin)
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
            setenv('SUBJECTS_DIR', fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211', 'derivatives'))
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
            switch lower(ipr.physiology)
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
            theData = QuadraticAerobicGlycolysisKit.constructData( ...
                'subjectsExpr', ipr.subjectsExpr, ...
                'tracer', tracer, ...
                'metric', metric, ...
                'region', region); 
            these = cell(size(theData));
            if ipr.Nthreads > 1                
                parfor (p = 1:length(theData), ipr.Nthreads)
                    try
                        these{p} = construction(theData(p)); %#ok<PFBNS>
                    catch ME
                        handwarning(ME)
                    end
                end
            elseif ipr.Nthreads == 1
                for p = length(theData):-1:1
                    try
                        these{p} = construction(theData(p)); % RAM ~ 3.3 GB
                    catch ME
                        handwarning(ME)
                    end
                end
            end            
            popd(pwd1);
        end
        function this = constructCbfByRegion(varargin)
            %% CONSTRUCTCBFBYREGION
            %  @param required imagingData is mlpipeline.ImagingData.
            %  @return cbf on filesystem.
            
            this = mlvg.QuadraticAerobicGlycolysisKit(varargin{:});            
            Region = [upper(this.imagingData.regionTag(1)) this.imagingData.regionTag(2:end)];

            pwd0 = pushd(this.imagingData.subjectPath);
            fs_ = this.(['buildFsBy' Region])();             
            cbf_ = this.fs2cbf(fs_);
            cbf_.save() % save ImagingContext2            
            popd(pwd0);
        end
        function this = constructCbvByRegion(varargin)
            %% CONSTRUCTCBVBYREGION
            %  @param required imagingData is mlpipeline.ImagingData.
            %  @return cbv on filesystem.
            
            this = mlvg.QuadraticAerobicGlycolysisKit(varargin{:});
            Region = [upper(this.imagingData.regionTag(1)) this.imagingData.regionTag(2:end)];

            pwd0 = pushd(this.imagingData.subjectPath);            
            vs_ = this.(['buildVsBy' Region])(); 
            cbv_ = this.vs2cbv(vs_);
            cbv_.save() % save ImagingContext2
            popd(pwd0);
        end 
        function this = constructCmro2ByRegion(varargin)
            %% CONSTRUCTCMRO2BYREGION
            %  @param required imagingData is mlpipeline.ImagingData.
            %  @return cmro2 on filesystem.
            %  @return oef on filesystem.
            
            this = mlvg.QuadraticAerobicGlycolysisKit(varargin{:});            
%             this.constructPhysiologyDateOnly('cbf', ...
%                 'subjectFolder', this.imagingData.subjectFolder, ...
%                 'region', this.imagingData.regionTag, ...
%                 'imagingData', this.imagingData)
%             this.constructPhysiologyDateOnly('cbv', ...
%                 'subjectFolder', this.imagingData.subjectFolder, ...
%                 'region', this.imagingData.regionTag, ...
%                 'imagingData', this.imagingData)
            Region = [upper(this.imagingData.regionTag(1)) this.imagingData.regionTag(2:end)];

            pwd0 = pushd(this.imagingData.subjectPath);             
            os_ = this.(['buildOsBy' Region])();   
            cbf_ = this.cbfOnAtlas( ...
                'typ', 'mlfourd.ImagingContext2', ...
                'dateonly', false, ...
                'tags', [this.blurTag this.regionTag]);
            [cmro2_,oef_] = this.os2cmro2(os_, cbf_, this.model);
            cmro2_.save() % save ImagingContext2
            oef_.save()            
            popd(pwd0);
        end 
        function dat = constructData(varargin)
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'subjectsExpr', 'sub-*', @istext)
            addParameter(ip, 'sessionsExpr', 'ses-*/pet', @istext)
            addParameter(ip, 'tracer', '', @istext)
            addParameter(ip, 'metric', '', @istext)
            addParameter(ip, 'region', 'voxels', @istext)
            parse(ip, varargin{:})
            ipr = ip.Results;

            idx = 1;
            scans = glob(fullfile( ...
                getenv('SUBJECTS_DIR'), ipr.subjectsExpr, ipr.sessionsExpr, ...
                sprintf('*_trc-%s_proc-dyn_pet_on_T1w.nii.gz', ipr.tracer)))';
            for s = scans
                dat_ = mlvg.Ccir1211Mediator(s{1}); %#ok<AGROW>
                dat_.metric = ipr.metric;
                dat_.regionTag = ipr.region;
                dat(idx) = dat_; %#ok<AGROW> 
                idx = idx + 1;
            end            
        end 
        function ic = constructPhysiologyDateOnly(varargin)
            %% e.g., constructPhysiologyDateOnly('cbv', 'imagingData', obj) % pwd includes sub-108293            
            %  e.g., constructPhysiologyDateOnly('cbv', 'imagingData', obj, ...
            %                                    'workpath', '/path/to/sub-108293/resampling_restricted')
            %  e.g., constructPhysiologyDateOnly('cbv', 'imagingData', obj, ...
            %                                    'workpath', '/path/to/sub-108293/resampling_restricted', ...
            %                                    'filepatt', 'ocdt2021*_111_voxels.4dfp.hdr')
            
            import mlvg.QuadraticAerobicGlycolysisKit
            
            reg = mlvg.Ccir1211Registry.instance();
            ip = inputParser;
            addRequired(ip, 'physiology', @istext)
            addParameter(ip, 'imagingData', [], @(x) isa(x, 'mlpipeline.ImagingData'))
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
                if contains(g{ig}, ipr.imagingData.defects)
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
        function ps = petPointSpread()
            ps = 3.57 * sqrt(0.5);
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
    
    properties (Dependent)
        atlasTag
        blurTag
        scanFolder % e.g., resampling_restricted
        scanPath
        registry
        regionTag
        subjectFolder % e.g., sub-108293
        subjectPath
        tags
    end 

    methods % GET
        function g = get.atlasTag(~)
            g = mlvg.Ccir1211Registry.instance.atlasTag;
        end
        function g = get.blurTag(~)
            g = mlvg.Ccir1211Registry.instance.blurTag;
        end
        function g = get.scanFolder(this)
            g = this.imagingData.scanFolder;
        end  
        function g = get.scanPath(this)
            g = this.imagingData.scanPath;
        end  
        function g = get.registry(this)
            g = mlvg.CCIR1211Registry.instance();
        end
        function g = get.regionTag(this)
            g = this.imagingData.regionTag;
        end
        function g = get.subjectFolder(this)
            g = this.imagingData.subjectFolder;
        end
        function g = get.subjectPath(this)
            g = this.imagingData.subjectPath;
        end
        function g = get.tags(this)
            g = strcat(this.blurTag, this.regionTag);
        end
    end

    methods
        function fs_ = buildFsByVoxels(this, varargin)
            %% BUILDFSBYVOXELS
            %  @return fs in R^4 as mlfourd.ImagingContext2, without saving to filesystems.  
            
            import mloxygen.QuadraticNumericRaichle1983
            
            ensuredir(this.scanPath);
            pwd0 = pushd(this.scanPath);  
                                    
            icv = this.dlicv();
            devkit = mlpet.ScannerKit.createFromSession(this.imagingData);             
            scanner = devkit.buildScannerDevice();
            scannerBrain = scanner.volumeAveraged(icv);
            arterial = this.buildAif(devkit, scanner, scannerBrain);
            
            fs_ = copy(icv.imagingFormat);
            fs_.filepath = this.scanPath;
            ic = this.fsOnAtlas(tags=this.tags);
            fs_.fileprefix = ic.fileprefix;

            % solve Raichle
            fprintf('%s\n', datestr(now))
            fprintf('starting mlraichle.QuadraticAerobicGlycolysisKit.buildFsByVoxels\n')
            raichle = QuadraticNumericRaichle1983.createFromDeviceKit( ...
                devkit, ...
                'scanner', scanner, ...
                'arterial', arterial, ...
                'roi', icv);  
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
            
            ensuredir(this.scanPath);
            pwd0 = pushd(this.scanPath);  
                                    
            icv = this.dlicv();
            devkit = mlpet.ScannerKit.createFromSession(this.imagingData);            
            scanner = devkit.buildScannerDevice(); 
            scannerBrain = scanner.volumeAveraged(icv); 
            arterial = this.buildAif(devkit, scanner, scannerBrain);
            
            os_ = copy(icv.imagingFormat);
            os_.filepath = this.scanPath;
            ic = this.osOnAtlas(tags=this.tags);
            os_.fileprefix = ic.fileprefix;

            % solve Mintun
            fprintf('%s\n', datestr(now))
            fprintf('starting mlraichle.QuadraticAerobicGlycolysisKit.buildOsByVoxels\n')
            mintun = QuadraticNumericMintun1984.createFromDeviceKit( ...
                devkit, ...
                'scanner', scanner, ...
                'arterial', arterial, ...
                'roi', icv);  
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
            
            ensuredir(this.scanPath);
            pwd0 = pushd(this.scanPath);                                    
            
            icv = this.dlicv();
            devkit = mlpet.ScannerKit.createFromSession(this.imagingData);             
            scanner = devkit.buildScannerDevice();
            scannerBrain = scanner.volumeAveraged(icv);
            arterial = this.buildAif(devkit, scanner, scannerBrain);
            
            vs_ = copy(icv.imagingFormat);
            vs_.filepath = this.scanPath;
            obj = this.vsOnAtlas(tags=this.tags);
            vs_.fileprefix = obj.fileprefix;

            % solve Martin
            fprintf('%s\n', datestr(now))
            fprintf('starting mlraiche.QuadraticAerobicGlycolysisKit.buildVsByVoxels\n')
            martin = QuadraticNumericMartin1987.createFromDeviceKit( ...
                devkit, ...
                'scanner', scanner, ...
                'arterial', arterial, ...
                'roi', icv);
            martin = martin.solve();  
            this.model = martin;

            % insert Martin solutions into fs
            vs_.img = martin.vs('typ', 'single');

            vs_ = mlfourd.ImagingContext2(vs_);
            popd(pwd0);
        end
        function obj = dlicv(this, varargin)
            obj = this.imagingData.bids.dlicv_ic;
        end
        function obj = metricOnAtlas(this, metric, varargin)
            %% METRICONATLAS appends fileprefixes with information from this.dataAugmentation
            %  @param required metric is char.
            %  @param datetime is datetime or char, .e.g., '20200101000000' | ''.
            %  @param dateonly is logical.
            %  @param tags is char, e.g., 'b43_wmparc1', default ''.
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'metric', @ischar)
            addParameter(ip, 'datetime', this.imagingData.datetime, @(x) isdatetime(x) || ischar(x))
            addParameter(ip, 'dateonly', false, @islogical)
            addParameter(ip, 'tags', '', @ischar)
            parse(ip, metric, varargin{:})
            ipr = ip.Results;

            try
                g = glob(fullfile(this.imagingData.scanPath, ...
                    sprintf('*_%s_*%s*.nii.gz', metric, ipr.tags)));
                if ~isempty(g)
                    obj = mlfourd.ImagingContext2(g{1});
                    return
                end
            catch ME
                handexcept(ME)
            end
            
            if ~isempty(ipr.tags)
                ipr.tags = strip(ipr.tags, "_");
            end   
            if ischar(ipr.datetime)
                adatestr = ipr.datetime;
            end
            if isdatetime(ipr.datetime)
                if ipr.dateonly
                    adatestr = ['ses-' datestr(ipr.datetime, 'yyyymmdd') '000000'];
                else
                    adatestr = ['ses-' datestr(ipr.datetime, 'yyyymmddHHMMSS')];
                end
            end
            
            s = this.imagingData.bids.filename2struct(this.imagingData.imagingContext.fqfn);
            s.ses = adatestr;
            s.modal = ipr.metric;
            s.tag = ipr.tags;
            fqfn = this.imagingData.bids.struct2filename(s);
            obj = mlfourd.ImagingContext2(fqfn);
        end	  
    end

    %% PROTECTED

    methods (Access = {?mlpet.AbstractAerobicGlycolysisKit2, ?mlvg.QuadraticAerobicGlycolysisKit})
        function this = QuadraticAerobicGlycolysisKit(varargin)
            %% QUADRATICAEROBICGLYCOLYSISKIT 
            %  Args:
            %      imagingData (ImagingData): session-specific objects
            %      aifMethods: containers.Map
            
            this = this@mlpet.AbstractAerobicGlycolysisKit2(varargin{:})
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
