classdef Ccir1211Bids < handle & mlsiemens.VisionBids
    %% Key functions:
    %  unpack_cnda() uses dcm2fileprefix(), info2*(), save_json().
    %  
    %  Created 23-Jan-2022 22:39:48 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.11.0.1837725 (R2021b) Update 2 for MACI64.  Copyright 2022 John J. Lee.
    
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
        function this = create(varargin)
            this = mlvg.Ccir1211Bids(varargin{:});
        end
        function [fn,info] = dcm2fileprefix(dcm)
            %  Params:
            %      dcm (filename): understood by dicominfo().
            %  Returns:
            %      fn (filename): for use by dcm2niix().
            %      info (struct): from dicominfo().

            assert(~isempty(dcm))
            if iscell(dcm)
                dcm = dcm{1};
            end
            info = dicominfo(dcm);
            sid = info.PatientName.FamilyName;
            dt = strcat(info.AcquisitionDate, strtok(info.AcquisitionTime, '.'));
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
        end
        function pth = info2destination(pth, info)
            %  Params:
            %      pth (folder): e.g., /Users/jjlee/Singularity/CCIR_02111/sourcedata/sub-108293/{anat,}
            %      info (struct): from dicominfo().
            %  Returns:
            %      pth: e.g., /Users/jjlee/Singularity/CCIR_02111/sourcedata/sub-108293/{anat,fmap,func,pet}

            assert(contains(pth, 'sub-'))
            [pth1,fold] = fileparts(pth);            
            if contains(fold, mlvg.Ccir1211Bids.BIDS_MODALITIES)
                pth = pth1;
            end
            if contains(info.Modality, 'MR')
                if contains(lower(info.SeriesDescription), 'fieldmap')
                    pth = fullfile(pth, 'fmap');
                    ensuredir(pth)
                    return
                end
                if contains(lower(info.SeriesDescription), 'rest') || ...
                        contains(lower(info.SeriesDescription), 'ase')
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
            %      destination (folder):  e.g., /Users/jjlee/Singularity/CCIR_02111/sourcedata/sub-108293/anat

            import mlvg.Ccir1211Bids.*

            ip = inputParser;
            addParameter(ip, 'download_folder', pwd, @isfolder)
            addParameter(ip, 'destination', pwd, @isfolder)
            parse(ip, varargin{:})
            ipr = ip.Results;

            for series = globFoldersT(fullfile(ipr.download_folder, '*'))
                try
                    dcms = glob(fullfile(series{1}, 'DICOM', '*.dcm'));
                    [fileprefix,info] = dcm2fileprefix(dcms);
                    if contains(lower(info.SeriesDescription), 'setter')
                        continue
                    end
                    dest = info2destination(ipr.destination, info);
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
    end

    properties (Constant)
        BIDS_MODALITIES = {'anat' 'fmap' 'func' 'mri' 'pet'}
        PROJECT_FOLDER = 'CCIR_01211'
        SURFER_VERSION = '7.2.0'
    end

    properties
        flair_toglob
        pet_dyn_toglob
        pet_static_toglob
        t1w_toglob
        t2w_toglob
        tof_toglob
    end

    methods
        function j = json(this)
            j = this.json_;
        end
        function r = registry(this)
            r = this.registry_;
        end

        function this = Ccir1211Bids(varargin)
            %  Args:
            %      destinationPath (folder): will receive outputs.  Specify project ID & subject ID.
            %      projectPath (folder): belongs to a CCIR project.  
            %      subjectFolder (text): is the BIDS-adherent string for subject identity.
            %      subjectFolder (text): is the BIDS-adherent string for subject identity.
            
            this = this@mlsiemens.VisionBids(varargin{:})          

            this.flair_toglob = fullfile(this.sourceAnatPath, 'sub-*_3D_FLAIR_Sag.nii.gz');
            this.pet_dyn_toglob = fullfile(this.sourcePetPath, 'sub-*_proc-dyn_pet.nii.gz');
            this.pet_static_toglob = fullfile(this.sourcePetPath, 'sub-*_proc-static_pet.nii.gz');
            this.t1w_toglob = fullfile(this.sourceAnatPath, 'sub-*_T1w_MPR_vNav_4e_RMS.nii.gz');
            this.t2w_toglob = fullfile(this.sourceAnatPath, 'sub-*_T2w_SPC_vNava.nii.gz');
            this.tof_toglob = fullfile(this.sourceAnatPath, 'sub-*_tof_fl3d_tra_p2_multi-slab.nii.gz');

            this.json_ = mlvg.Ccir1211Json();
            this.registry_ = mlvg.Ccir1211Registry.instance();
        end
    end    

    %% PROTECTED

    properties (Access = protected)
        json_
        registry_
    end

    methods (Access = protected)
        function that = copyElement(this)
            that = copyElement@matlab.mixin.Copyable(this);
        end
    end

    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
