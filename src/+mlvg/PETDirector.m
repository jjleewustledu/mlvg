classdef PETDirector < handle
    %% directs building PET-relevant representations, 
    %  progressively updating internal builders
    %  
    %  Created 05-Sep-2022 15:50:41 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.12.0.2039608 (R2022a) Update 5 for MACI64.  Copyright 2022 John J. Lee.

    properties (Dependent)
        product
        study_builder
        reconstruction_builder
        workdir
    end

    methods

        %% GET

        function g = get.product(this)
            g = [];
        end
        function g = get.study_builder(this)
            g = this.study_builder_;
        end
        function g = get.reconstruction_builder(this)
            g = this.reconstruction_builder_;
        end
        function g = get.workdir(this)
            g = this.workdir_;
        end

        %%

        function this = build_filesystem(this, varargin)
            this.study_builder_.build_filesystem(varargin{:});
        end
        function this = build_unpacked(this, varargin)
            this.study_builder_.build_unpacked(varargin{:});
        end
        function this = build_reconstruction(this, varargin)
            call(this.reconstruction_builder_, varargin{:});
        end
        function this = build_bids(this, varargin)
            %% $SINGULARITY_HOME/CCIR_01211/
            %  {rawdata,sourcedata,derivatives}/sub-<id>/ses-<date>/
            %  {anat,fmap,func,pet}/*.nii.gz

            call(this.bids_builder_, varargin{:});
        end
        function this = build_calibration(this, varargin)
        end
        function this = build_imaging(this, varargin)
        end
        function this = build_input_function(this, varargin)
        end
        function this = build_model(this, varargin)
        end
        function this = build_quality_assurance(this, varargin)
        end
        function this = build_reports(this, varargin)
        end
        function this = build_tracer(this, varargin)
        end
        function this = call(this, varargin)
            %% CALL typical constructions, representations

            this.build_filesystem(varargin{:});

            % rawdata:
            % -> build tracer & imaging representations
            this.build_unpacked(varargin{:});
            this.build_tracer(varargin{:});
            this.build_reconstruction(varargin{:});

            % sourcedata:  
            % -> generate BIDS
            this.build_bids(varargin{:});

            % derivatives:
            % -> build device & data representations
            this.build_imaging(varargin{:});
            this.build_input_function(varargin{:});

            % -> cross-calibrate devices, align configurations & spacetimes
            this.build_calibration(varargin{:});

            % -> build scientific models
            this.build_model(varargin{:});
            this.build_reports(varargin{:});
        end
        function self_test(this)
            %error('mlvg:NotImplementedError', 'PETDirector.self_test')
        end

        function this = PETDirector(varargin)
            %% PETDIRECTOR 
            %  Args:
            %      workdir (text): e.g., "rawdata/sub-108287", 
            %                            "sourcedata/sub-108287",
            %                            "derivatives/sub-108*"
            %                            fullfile(getenv("SINGULARITY_HOME"), ...
            %                                    "Singularity/CCIR_01211/derivatives/sub-108287").
            %      reconstruction_builder (text):  e.g., "niftypet", "jsrecon"

            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, "workdir", pwd, @isfolder);
            addParameter(ip, "reconstruction_builder", "", @istext)
            parse(ip, varargin{:});
            ipr = ip.Results;
            this.workdir_ = ipr.workdir;
            this.reconstruction_builder_ = ipr.reconstruction_builder;

            % builder design patterns usefully retain reference to director
            this.select_study_builder(varargin{:});
            this.select_reconstruction_builder(varargin{:});
        end
    end

    %% PROTECTED

    properties % (Access = protected)
        study_builder_
        bids_builder_
        reconstruction_builder_
        imaging_builder_
        tracer_builder_
        input_function_builder_
        calibration_builder_
        model_builder_

        workdir_
    end

    methods (Access = protected)
        function this = select_study_builder(this, varargin)
            if contains(this.workdir, 'CCIR_00754') || contains(this.workdir, 'CCIR_00559')
                this.study_builder_ = mlraichle.Ccir559754Builder(this, varargin{:});
                return
            end
            if contains(this.workdir, 'CCIR_00993')
                this.study_builder_ = mlan.Ccir993Builder(this, varargin{:});
                return
            end
            if contains(this.workdir, 'CCIR_01211')
                this.study_builder_ = mlvg.Ccir1211Builder(this, varargin{:});
                return
            end

            error('mlvg:ValueError', '%s: %s not supported', ...
                clientname(), this.workdir)
        end
        function this = select_reconstruction_builder(this, varargin)
            if contains(this.reconstruction_builder_, 'niftypet', 'IgnoreCase', true)
                this.reconstruction_builder_ = mlnipet.NiftypetBuilder(this, varargin{:});
                return
            end
            if contains(this.reconstruction_builder_, 'jsrecon', 'IgnoreCase', true)
                this.reconstruction_builder_ = mlsiemens.JSReconBuilder(this, varargin{:});
                return
            end

            if contains(this.workdir, 'CCIR_00754') || contains(this.workdir, 'CCIR_00559')
                this.reconstruction_builder_ = mlnipet.NiftypetBuilder(this, varargin{:});
                return
            end
            if contains(this.workdir, 'CCIR_00993')
                this.reconstruction_builder_ = mlnipet.NiftypetBuilder(this, varargin{:});
                return
            end
            if contains(this.workdir, 'CCIR_01211')
                this.reconstruction_builder_ = mlsiemens.JSReconBuilder(this, varargin{:});
                return
            end

            error('mlvg:ValueError', '%s: %s & %s not supported', ...
                clientname(), this.workdir, class(this.reconstruction_builder_))
        end
        function this = select_bids_builder(this, varargin)
        end
        function this = select_imaging_builder(this, varargin)
        end
        function this = select_tracer_builder(this, varargin)
        end
        function this = select_input_function_builder(this, varargin)
        end
        function this = select_calibration_builder(this, varargin)
        end
        function this = select_model_builder(this, varargin)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
