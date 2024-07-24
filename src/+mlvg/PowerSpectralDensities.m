classdef PowerSpectralDensities < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 10-Jul-2024 12:35:46 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 24.1.0.2628055 (R2024a) Update 4 for MACA64.  Copyright 2024 John J. Lee.
    

    properties
        pet
        t0
    end

    properties (Dependent)
        adj_times
        img
        sampling_freq  % Hz

        iFV  % tentatively just entire 4th ventricle
        brainstem
        thalamus
        white
        white_precuneus
        grey
        grey_precuneus        
    end

    methods %% GET
        function g = get.adj_times(this)
            t1 = this.t0 + 1;
            g = this.pet.json_metadata.times(t1:end);            
        end

        function g = get.img(this)
            if ~isempty(this.img_)
                g = this.img_;
                return
            end

            t1 = this.t0 + 1;
            this.img_ = this.pet.imagingFormat.img(:,t1:end);
            g = this.img_;
        end

        function g = get.sampling_freq(this)
            dt = this.adj_times(2) - this.adj_times(1);
            g = 1/dt;
        end

        function g = get.iFV(this)
            g = this.img(11, :)';
        end
        function g = get.brainstem(this)
            g = mean(this.img([12, 17, 32], :), 1)';
        end
        function g = get.thalamus(this)
            g = mean(this.img([6, 25], :), 1)';
        end
        function g = get.white(this)
            g = mean(this.img(42:109, :), 1)';
        end
        function g = get.white_precuneus(this)
            g = mean(this.img([65, 99], :), 1)';
        end
        function g = get.grey(this)
            g = mean(this.img(110:309, :), 1)';
        end
        function g = get.grey_precuneus(this)
            g = mean(this.img([204:208, 307:309], :), 1)';
        end
    end

    methods
        function this = PowerSpectralDensities(pet, opts)
            arguments
                pet {mustBeNonempty}
                opts.t0 double = 0
            end

            this.pet = mlfourd.ImagingContext2(pet);
            assert(contains(this.pet.fileprefix, "ParcSchaeffer"))
            this.t0 = opts.t0;
        end
    end

    %% PRIVATE

    properties (Access=private)
        img_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
