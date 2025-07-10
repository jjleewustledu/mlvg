classdef Ccir1211Bids < handle & mlsiemens.BiographBids
    %% Key functions:
    %  unpack_cnda() uses dcm2fileprefix(), info2*(), save_json().
    %  
    %  Created 23-Jan-2022 22:39:48 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.11.0.1837725 (R2021b) Update 2 for MACI64.  Copyright 2022 John J. Lee.
    
    properties (Constant)
        BIDS_MODALITIES = {'anat' 'fmap' 'func' 'mri' 'pet'}
        DLICV_TAG = 'DLICV'
        PROJECT_FOLDER = 'CCIR_01211'
        SURFER_VERSION = '7.3.2'
    end

    methods
        function this = Ccir1211Bids(varargin)
            %  Args:
            %      destinationPath (folder): will receive outputs.  Specify project ID & subject ID.
            %      projectPath (folder): belongs to a CCIR project.  
            %      subjectFolder (text): is the BIDS-adherent string for subject identity.
            %      subjectFolder (text): is the BIDS-adherent string for subject identity.

            this = this@mlsiemens.BiographBids(varargin{:})
            try
                this.json_ = mlvg.Ccir1211Json();
            catch ME
                fprintf("%s: safely ignoring mlvg.Ccir1211Json\n", stackstr());
            end
        end
        function r = registry(~)
            r = mlvg.Ccir1211Registry.instance();
        end
    end    

    methods (Static)
        function bids_to_4dfp(varargin)
            ip = inputParser;
            addRequired(ip, 'subjectPath', pwd, @isfolder)
            parse(ip, varargin{:})
            ipr = ip.Results;

            % ensure subjectPath
            c = fileparts2cell(ipr.subjectPath);
            [~,pos] = max(contains(c, 'sub-'));
            ipr.subjectPath = fullfile(filesep, c{1:pos});

            % ensure dataPath ~ resampling_restricted
            dataPath = fullfile(ipr.subjectPath, 'resampling_restricted', '');
            ensuredir(dataPath)

            % generate 4dfp
            for nii = globFoldersT(fullfile(ipr.subjectPath, 'pet', 'ses-*', 'sub-*_ses-*_proc-dyn_pet_on_T1w.nii.gz'))
                ic = mlfourd.ImagingContext2(nii{1});
                ic.selectFourdfpTool();
                ic.filepath = dataPath;
                [~,fp] = myfileparts(nii{1});
                re = regexp(fp, 'sub-\d{6}_ses-(?<dt>\d{14})_trc-(?<trc>\w+)_proc-dyn_pet_on_T1w', 'names');
                ic.fileprefix = strcat(re.trc, 'dt', re.dt);
                ic.save();
                ic.timeAveraged();
                ic.save();
            end
        end
        
        function convert_metcalf_niftis(petdir, destdir, opts)
            %% identifies console recons, copying to destination with
            %  renaming .nii.gz and .json, and adding json fields "taus", "timesMid".

            arguments
                petdir {mustBeFolder} = pwd
                destdir {mustBeFolder} = pwd
                opts.tracer_patt {mustBeTextScalar} = "FDG Dynamic"  % in json
                opts.tracer {mustBeTextScalar} = "fdg"
                opts.tag {mustBeTextScalar} = "consoleDynamic"
            end

            pwd0 = pushd(petdir);
            try

                % inspect json files to match tracer_patt

                jfiles = [];
                for jfile = mglob(fullfile(petdir, "*.json"))
                    j = jsondecode(fileread(jfile));
                    if strcmp(j.SeriesDescription, opts.tracer_patt)
                        jfiles = [jfiles, jfile]; %#ok<AGROW>
                    end
                end

                % generate new filenames

                jfiles1 = [];
                for jfile = jfiles
                    re = regexp(jfile, "(?<sub>sub-\d{6})/(?<ses>ses-\d{8})", "names");
                    sub = re.sub;
                    ses = re.ses;
                    fqfp = fullfile( ...
                        destdir, sprintf("%s_%s_trc-%s_proc-%s", sub, ses, opts.tracer, opts.tag));
                    jfiles1 = [jfiles1, fqfp + ".json"]; %#ok<AGROW>
                end

                % copy json and nii.gz to destination

                nfiles = strrep(jfiles, ".json", ".nii.gz");
                nfiles1 = strrep(jfiles1, ".json", ".nii.gz");
                for jidx = 1:length(jfiles)
                    copyfile(jfiles(jidx), jfiles1(jidx));
                    copyfile(nfiles(jidx), nfiles1(jidx));
                end

                % for FDG Dynamic add taus, timesMid to json, updating json on filesystem

                if contains(opts.tracer_patt, "Dynamic")
                    for jidx = 1:length(jfiles1)
                        jfile1 = jfiles1(jidx);
                        j1 = jsondecode(fileread(jfile1));
                        j1.taus = j1.FrameDuration;
                        j1.timesMid = j1.FrameReferenceTime;
                        writelines(jsonencode(j1, PrettyPrint=true), jfile1);
                    end
                end

            catch ME
                handwarning(ME)
            end
            popd(pwd0);
        end

        function this = create(varargin)
            this = mlvg.Ccir1211Bids(varargin{:});
        end
        function [fn,info,ses_fold] = dcm2fileprefix(dcm)
            %  Params:
            %      dcm (filename): understood by dicominfo().
            %  Returns:
            %      fn (filename): for use by dcm2niix().
            %      info (struct): from dicominfo().
            %      ses_fold (text):  e.g., "ses-20250101".

            assert(~isempty(dcm))
            if iscell(dcm)
                dcm = dcm{1};
            end
            info = dicominfo(dcm);
            sid = info.PatientName.FamilyName;
            if isfield(info, 'AcquisitionDateTime')
                dt = strtok(info.AcquisitionDateTime, '.');
            elseif isfield(info, 'AcquisitionDate') && isfield(info, 'AcquisitionTime')
                dt = strcat(info.AcquisitionDate, strtok(info.AcquisitionTime, '.'));
            else
                error('mlvg:RuntimeError', stackstr())
            end
            switch info.Modality
                case 'CT'
                    fn = sprintf('sub-%s_ses-%s_kvp-%g_ct', sid, dt, info.KVP);
                case 'PT'
                    [trc,proc] = mlvg.Ccir1211Bids.info2pet_fields(info);
                    fn = sprintf('sub-%s_ses-%s_trc-%s_proc-%s_pet', sid, dt, trc, proc);
                case 'MR'
                    pulses = strrep(info.SeriesDescription, ' ', '_');
                    fn = sprintf('sub-%s_ses-%s_%s', sid, dt, pulses);
                otherwise 
                    dsc = strrep(info.SeriesDescription, ' ', '_');
                    fn = sprintf('sub-%s_ses-%s_dsc-%s-%s_%s', sid, dt, dsc, info.Modality);
            end

            ses_fold = ['ses-', dt(1:8)];  % yyyyMMdd
        end
        function pth = info2destination(pth, info)
            %  Params:
            %      pth (folder): e.g., /Users/jjlee/Singularity/CCIR_02111/sourcedata/sub-108293/
            %      info (struct): from dicominfo().
            %  Returns:
            %      pth: e.g., /Users/jjlee/Singularity/CCIR_02111/sourcedata/sub-108293/ses-yyyyMMdd/{anat,fmap,func,pet}/

            assert(contains(pth, 'sub-'))
            [pth1,fold] = fileparts(pth);            
            if any(contains(fold, mlvg.Ccir1211Bids.BIDS_MODALITIES))
                pth = pth1;
            end
            if contains(info.Modality, 'MR')
                if contains(lower(info.SeriesDescription), 'localizer')
                    pth = fullfile(pth, 'localizer');
                    ensuredir(pth)
                    return
                end
                if contains(lower(info.SeriesDescription), 'fieldmap')
                    pth = fullfile(pth, 'fmap');
                    ensuredir(pth)
                    return
                end
                if contains(lower(info.SeriesDescription), 'dti')
                    pth = fullfile(pth, 'dwi');
                    ensuredir(pth)
                    return
                end
                if contains(lower(info.SeriesDescription), 'rest') || ...
                        contains(lower(info.SeriesDescription), 'ase') || ...
                        contains(lower(info.SeriesDescription), 'qsm')
                    pth = fullfile(pth, 'func');
                    ensuredir(pth)
                    return
                end
                pth = fullfile(pth, 'anat');
                ensuredir(pth)
                return
            end
            if contains(info.Modality, 'CT') || contains(info.Modality, 'PT')
                pth = fullfile(pth, 'pet');
                ensuredir(pth)
            end
            % pth == /path/to/sub-123456 for unknown dicoms
        end
        function txt = info2json(info)
            for f = asrow(fields(info))
                if numel(info.(f{1})) > 100
                    large = info.(f{1}); 
                    info.(f{1}) = large(1:20);
                end
            end
            txt = jsonencode(info, 'PrettyPrint', true);
            txt = strrep(txt, '\', '_');
        end
        function [trc,proc] = info2pet_fields(info)
            assert(isstruct(info))
            dsc = lower(info.SeriesDescription);
            trc = 'unknown';
            if contains(dsc, 'co') || contains(dsc, 'oc')
                trc = 'oc';
            end
            if contains(dsc, 'oo') || contains(dsc, 'oxygen')
                trc = 'oo';
            end
            if contains(dsc, 'ho') || contains(dsc, 'water')
                trc = 'ho';
            end
            if contains(dsc, 'fdg')
                trc = 'fdg';
            end
            if strcmp(trc, 'unknown')
                try
                    radio = lower(info.RadiopharmaceuticalInformationSequence.Item_1.Radiopharmaceutical);
                    if contains(radio, 'co') || contains(radio, 'oc')
                        trc = 'oc';
                    end
                    if contains(radio, 'oo') || contains(radio, 'oxygen')
                        trc = 'oo';
                    end
                    if contains(radio, 'ho') || contains(radio, 'water')
                        trc = 'ho';
                    end
                    if contains(radio, 'fdg')
                        trc = 'fdg';
                    end
                    if strcmp(trc, 'unknown')
                        trc = lower(radio);
                    end
                catch ME
                    handwarning(ME)
                end
            end
            proc = 'unknown';            
            if contains(dsc, 'static')
                proc = 'static';
            end
            if contains(dsc, 'dynamic')
                proc = 'dyn';
            end
            if contains(dsc, 'nac')
                proc = 'nac';
            end
        end
        function tf = isdynamic(obj)
            ic = mlfourd.ImagingContext2(obj);
            re = regexp(ic.fileprefix, '\S+_proc-(?<dyn>\w+)_\S+', 'names');
            tf = contains(re.dyn, 'dyn');
        end
        function tf = isnac(obj)
            ic = mlfourd.ImagingContext2(obj);
            re = regexp(ic.fileprefix, '\S+_proc-(?<nac>\w+)_\S+', 'names');
            tf = contains(re.nac, 'nac');
        end
        function tf = isstatic(obj)
            ic = mlfourd.ImagingContext2(obj);
            re = regexp(ic.fileprefix, '\S+_proc-(?<stat>\w+)_\S+', 'names');
            tf = contains(re.stat, 'stat');
        end
        function tr = obj2tracer(obj)
            ic = mlfourd.ImagingContext2(obj);
            re = regexp(ic.fileprefix, '\S+_trc-(?<tr>\w+)_\S+', 'names');
            tr = upper(re.tr);
        end
        function save_json(filename, j)
            if isstruct(j)
                j = mlvg.Ccir1211Bids.info2json(j);
            end
            fid = fopen(filename, 'w');
            fprintf(fid, j);
            fclose(fid); 
        end
        function unpack_cnda(varargin)
            %  Params:
            %      download_folder (folder):  e.g., /Users/jjlee/Downloads/CNDA_XNAT_Desktop_Client/cnda.wustl.edu/108293_MAG_20210421
            %      destination (folder):  e.g., /Users/jjlee/Singularity/CCIR_02111/sourcedata/sub-108293

            import mlvg.Ccir1211Bids.*

            ip = inputParser;
            addParameter(ip, 'download_folder', pwd, @isfolder)
            addParameter(ip, 'destination', pwd, @isfolder)
            addParameter(ip, 'modality_folder', '', @istext)
            parse(ip, varargin{:})
            ipr = ip.Results;

            for series = globFoldersT(fullfile(ipr.download_folder, '*'))
                try
                    dcms = glob(fullfile(series{1}, 'DICOM', '*.dcm'));
                    [fileprefix,info,ses_fold] = dcm2fileprefix(dcms);
                    if contains(lower(info.SeriesDescription), 'setter')
                        continue
                    end
                    dest = fullfile(ipr.destination, ses_fold);
                    dest = info2destination(dest, info);
                    if ~isemptytext(ipr.modality_folder)
                        if ~endsWith(dest, ipr.modality_folder)
                            continue
                        end
                    end
                    mlvg.Ccir1211Bids.dcm2niix(series{1}, ...
                        'f', fileprefix, ...
                        'o', dest, ...
                        'fourdfp', false);
                    save_json(fullfile(dest, strcat(fileprefix, '.dcm.json')), info);
                catch ME
                    handwarning(ME)
                end
            end
        end
        function unpack_metcalf_dicoms(dicoms_best_per_sub)
            arguments
                dicoms_best_per_sub {mustBeFile} = ...
                    "~/mnt/CHPC_scratch/Singularity/CCIR_01211/derivatives/dicoms_best_per_sub.mat"
            end

            import mlvg.Ccir1211Bids.*

            ld = load(dicoms_best_per_sub);
            for d = ld.dicoms1  
                % e.g.: /data/nil-bluearc/vlassenko/RAW_IMAGES/MRI/108007/MR2_20201203/DICOMS

                sub = "sub-" + find_108_substring(d);
                destination = fullfile("~/mnt/CHPC_scratch/Singularity/CCIR_01211/sourcedata", sub);
                ensuredir(destination);
                unpack_cnda( ...
                    download_folder=d, ...
                    destination=destination, ...
                    modality_folder="anat");
            end
        end
        function unpack_RAW_IMAGES()

            MIR = "/data/nil-bluearc/vlassenko/RAW_IMAGES/MRI";
            sourcedata = "/data/nil-bluearc/vlassenko/jjlee/Singularity/CCIR_01211/sourcedata";
            source_folders = readlines(fullfile(sourcedata, "source_folders.log"));
            subids = readlines(fullfile(sourcedata, "subids.log"));
        end

        %% HELPERS

        function digit_string = find_108_substring(unix_path)
            % Find all occurrences of /108xxx/ pattern
            pattern = '/108\d{3}/';
            matches = regexp(unix_path, pattern, 'match');

            if ~isempty(matches)
                % Extract the digits from the first match
                digit_string = extractBetween(matches{1}, '/', '/');
                digit_string = digit_string{1};  % extractBetween returns cell array
            else
                digit_string = '';
                warning('No 6-digit substring beginning with "108" found between fileseps');
            end
        end
    end

    %% PROTECTED

    methods (Access = protected)
        function that = copyElement(this)
            that = copyElement@matlab.mixin.Copyable(this);
        end
    end

    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
