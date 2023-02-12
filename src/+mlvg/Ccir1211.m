classdef Ccir1211 < handle
    %% line1
    %  line2
    %  
    %  Created 09-Jun-2022 15:34:35 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.12.0.1956245 (R2022a) Update 2 for MACI64.  Copyright 2022 John J. Lee.
    
    properties (Dependent)
        bids
        study
    end

    methods

        %% GET

        function g = get.bids(this)
            g = this.bids_;
        end
        function g = get.study(this)
            g = this.study_;
        end

        %%

        function this = Ccir1211(varargin)
            %% CCIR1211 
            %  Args:
            %      arg1 (its_class): Description of arg1.
            
            ip = inputParser;
            addParameter(ip, "arg1", [], @(x) true)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.bids_ = mlvg.Ccir1211Bids.create();
            this.study_ = mlvg.Ccir1211Registry.instance();
        end

        function call(this)
        end
        function call_sub(this, sub)
        end
        function mdl_agi = call_ses(this, ses)
            oxy = mlkinetics.OxygenMetabKit.instance();                
            tracers = this.make_tracer_iterator(oxy);
            mdl_oxy = {};
            m = 1;
            while tracers.hasNext()
                tra = tracers.next();
                sca = oxy.make_scanner(this, ses, tra);
                aif = oxy.make_input_function(this, sca);
                mdl_oxy{m} = oxy.make_model(this, sca, aif); %#ok<AGROW> 
                m = m + 1;
            end

            glc = mlkinetics.GlucoseMetabKit.instance();
            g = glob('*_trc-fdg_proc-dyn_pet.nii.gz');
            ic2 = mlfourd.ImagingContext2(g{1});
            tra = glc.make_tracer(this, ic2);
            sca = glc.make_scanner(this, ses, tra);
            aif = glc.make_input_function(this, sca);
            mdl_glc = glc.make_model(this, sca, aif);

            agi = mlkinetics.AerobicGlycolysisKit.instance();
            mdl_agi = agi.make_model(this, mdl_oxy, mdl_glc);
        end
        function iter = make_sub_iterator(this)
        end
        function iter = make_ses_iterator(this)
        end
        function iter = make_tracer_iterator(this, kinetics_kit)
            pwd0 = pushd(this.bids_.derivPetPath);
            iter = mlkinetics.TracerIterator(kinetics_kit, ...
                'patt', '*_trc-{oc,co,oo,ho}_proc-dyn_pet.nii.gz', ...
                'taus', @this.study_.taus);
            popd(pwd0);
        end
    end

    %% PRIVATE

    properties (Access = private)
        bids_
        study_
        taus_map_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
