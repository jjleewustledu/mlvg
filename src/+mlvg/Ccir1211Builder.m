classdef Ccir1211Builder < handle & mlpipeline.StudyBuilder
    %% Configures CCIR project 01211 for mlvg.PETDirector
    %  
    %  Created 06-Sep-2022 16:26:49 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.12.0.2039608 (R2022a) Update 5 for MACI64.  Copyright 2022 John J. Lee.
    
    properties (Dependent)
        % BIDS fully-qualified directories
        derivativesDir
        rawdataDir
        sourcedataDir
        studyDir
        studyFolder

        dicomExt
        studyJson
        workdir
    end

    methods % GET
        function g = get.derivativesDir(this)
            g = fullfile(this.studyDir, 'derivatives');
        end
        function g = get.dicomExt(~)
            g = '.dcm';
        end
        function g = get.rawdataDir(this)
            g = fullfile(this.studyDir, 'rawdata');
        end
        function g = get.sourcedataDir(this)
            g = fullfile(this.studyDir, 'sourcedata');
        end        
        function g = get.studyDir(this)
            g = fullfile(getenv('SINGULARITY_HOME'), this.studyFolder);
        end
        function g = get.studyFolder(~)
            g = 'CCIR_01211';
        end
        function g = get.studyJson(this)
            fn = fullfile(this.studyDir, 'constructed_20210225.json');
            g = jsondecode(fileread(fn));
        end
        function g = get.workdir(this)
            g = this.director_.workdir;
        end
    end
        
    methods        
        function this = build_bids(this, varargin)
            %% Builds filesystem with, e.g., $SINGULARITY_HOME/CCIR_*/
            %  {rawdata,sourcedata,derivatives}/sub-<id>/ses-<date>/
            %  {anat,fmap,func,pet}/*.nii.gz
        end
        function this = build_quality_assurance(this, varargin)
        end
        function this = build_unpacked(this)
            %% Unpacks *.zip found in workdir at or deeper than CCIR_01211/rawdata/sub-*/;
            %  creates *.nii.gz in CCIR_01211/sourcedata/sub-*/ses-*/modality.
            %  this.workdir must contain sub-*.

            assert(contains(this.workdir, 'sub-'), clientname(false, 2))

            % unzip *.zip
            zs = glob(fullfile(this.workdir, '**.zip'))';
            if ~isempty(zs)
                for z = zs

                    % pushd folder containing *.zip, e.g., sub-*/
                    pwd0 = pushd(fileparts(z{1}));
                    try
                        unzip(z{1});
                        %delete(z{1});
                    catch ME
                        handwarning(ME)
                    end
                    popd(pwd0);
                end
            end

            % dcm2niix from all DICOMs
            pwd0 = pushd(this.workdir);
            dcms = glob(fullfile(this.workdir, '**.dcm'))';

            % cell array := unique folders containing dcm
            dns = unique(cellfun(@(x) myfileparts(x), dcms, 'UniformOutput', false));
            dns = dns(~contains(dns, 'secondary'));
            subfold_ = this.build_subFolder(this.workdir);
            for d = dns
                try
                    dcm_ = glob(fullfile(d{1}, '**.dcm'));
                    sesfold_ = this.build_sesFolder(subfold_, dcm_{1});
                    modfold_ = this.build_modalityFolder(subfold_, sesfold_, dcm_{1});
                    targetSourceDir = fullfile(this.sourcedataDir, subfold_, sesfold_, modfold_);
                    ensuredir(targetSourceDir);

                    % *.dcm -> sub-*/ses-*/*.nii.gz
                    mlpipeline.Bids.dcm2niix(d{1}, 'o', targetSourceDir);

                    % clean-up targetSourceDir
                    pwd1 = pushd(targetSourceDir);
                    this.move_etc();
                    popd(pwd1);
                catch ME
                    handwarning(ME)
                end
            end
            popd(pwd0);
        end

        function this = Ccir1211Builder(director, varargin)
            %% CCIR1211BUILDER 
            %  Args:
            %      director mlvg.PETDirector : references related builders.

            arguments (Input)
                director mlvg.PETDirector
            end
            arguments (Input, Repeating)
                varargin % legacy support; https://www.mathworks.com/help/matlab/ref/varargin.html
            end

            this = this@mlpipeline.StudyBuilder(varargin{:});           
            
            this.director_ = director;
            this.registry_ = mlvg.Ccir1211Registry.instance();
        end
    end

    %% PROTECTED

    properties (Access = protected)
        director_
        registry_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
