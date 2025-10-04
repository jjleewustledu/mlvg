classdef T4Resolve < handle & mlsystem.IHandle
    %% run on machine.neuroimage.wustl.edu
    %  
    %  Created 04-Oct-2025 00:11:23 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.
    

    properties
        debug
    end

    properties (Dependent)
        pet
        t1w
        t1w_brain
        t1w_dlicv
        t1w_on_pet
    end

    methods  % get, set
        function g = get.pet(this)
            g = this.pet_;
            assert(isfile(g))
        end
        function g = get.t1w(this)
            assert(contains(this.t1w_.fileprefix, "orient-std"))
            g = this.t1w_;
            assert(isfile(g))
        end
        function g = get.t1w_brain(this)
            if ~isempty(this.t1w_brain_)
                g = this.t1w_brain_;
                return
            end

            % read from filesystem
            fqfn = this.t1w.fqfp + "_brain.nii.gz";
            if isfile(fqfn)
                this.t1w_brain_ = mlfourd.ImagingContext2(fqfn);
                g = this.t1w_brain_;
                return
            end

            % build and save nii.gz
            t1wb = this.t1w.imagingFormat;
            t1wb.img = single(t1wb.img) .* single(this.t1w_dlicv.imagingFormat.img);
            this.t1w_brain_ = mlfourd.ImagingContext2(t1wb);
            this.t1w_brain_.fqfn = fqfn;
            this.t1w_brain_.save();

            % build and save json
            in = this.t1w;            
            S = struct(stackstr(2), 'this.t1w_brain_ := this.t1w.imagingFormat.img .* this.t1w_dlicv.imagingFormat.img');
            this.jsonrecode(in, S, this.t1w_brain_);
            assert(isfile(this.t1w_brain_))
            g = this.t1w_brain_;
        end
        function g = get.t1w_dlicv(this)
            if ~isempty(this.t1w_dlicv_)
                g = this.t1w_dlicv_;
                return
            end

            fqfp = this.t1w.fqfp;
            this.t1w_dlicv_ = mlfourd.ImagingContext2(sprintf('%s_DLICV.nii.gz', fqfp));
            if ~isfile(this.t1w_dlicv_)
                mg = glob(fullfile(fileparts(fqfp), '*T1w*orient-std*_DLICV.nii.gz'));
                assert(~isempty(mg))
                mg = natsort(mg);
                this.t1w_dlicv_ = mlfourd.ImagingContext2(mg(end));
            end
            g = this.t1w_dlicv_;
        end
        function g = get.t1w_on_pet(this)
            filepath = strrep(this.pet.filepath, "sourcedata", "derivatives");
            fqfn = fullfile(filepath, "T1w_on_" + this.pet.filename);
            g = mlfourd.ImagingContext2(fqfn);
        end
    end

    methods
        function this = T4Resolve(opts)
            %% T4RESOLVE
            
            arguments
                opts.pet {mustBeNonempty}  % 3D representation, e.g. static, avgt, mipt
                opts.t1w {mustBeNonempty}
                opts.atl = fullfile(getenv("REFDIR"), "MNI152_T1_1mm.nii.gz")
                opts.blur double {mustBeScalarOrEmpty} = 3.5
                opts.debug logical = false
            end
            
            this.pet_ = mlfourd.ImagingContext2(opts.pet);
            this.t1w_ = opts.t1w;
            assert(3 == ndims(this.pet_))
            this.atl_ = mlfourd.ImagingContext2(opts.atl);
            this.blur_ = opts.blur;
            this.debug = opts.debug;
        end

        function this = flirt_t1w_to_pet(this, opts)
            arguments
                this mlvg.T4Resolve
                opts.noclobber logical = false
            end

            if opts.noclobber && isfile(this.t1w_on_pet)
                return
            end

            t1w__ = this.t1w_brain;
            flirt_pet = mlfsl.Flirt( ...
                'in', t1w__.fqfilename, ...
                'ref', this.pet.fqfilename, ...
                'out', this.niigz(this.t1w_on_pet), ...
                'omat', this.mat(this.t1w_on_pet), ...
                'bins', 4096, ...
                'cost', 'mutualinfo', ...
                'dof', 6, ...
                'searchrx', 45, ...
                'searchry', 45, ...
                'searchrz', 45, ...
                'interp', 'spline');
            flirt_pet.flirt();

            j0 = fileread(this.pet.fqfp + ".json");
            [~,j1] = flirt_pet.cost_final();
            jsonrecode(j0, j1, 'filenameNew', this.json(this.t1w_on_pet));
        end
        
        function this = resolve_t1w_to_pet(this, opts)
            arguments
                this mlvg.T4Resolve
                opts.noclobber logical = false
            end

            if opts.noclobber && isfile(this.t1w_on_pet)
                return
            end

            filepath = strrep(this.pet.filepath, "sourcedata", "derivatives");
            pwd0 = pushd(filepath);
            
            % resolve
            msks{1} = mlfourd.ImagingContext2('none.nii.gz');
            msks{2} = this.t1w_dlicv;
            imgs{1} = this.pet;
            imgs{2} = this.t1w_brain;
            t4rb = mlfourd.SimpleT4ResolveBuilder( ...
                'workpath', filepath, ...
                'maskForImages', msks, ...
                'theImages', imgs, ...
                'debug', this.debug);
            dispdbg(t4rb)
            t4rb = t4rb.resolve();
            
            % move files, write json
            movefile(this.niigz(t4rb.theImagesFinal{2}), this.niigz(this.t1w_on_pet));
            movefile(this.json(t4rb.theImagesFinal{2}), this.json(this.t1w_on_pet));    
            this.jsonrecode( ...
                this.t1w_on_pet, ...
                struct('image_activity', this.image_mass(this.t1w_on_pet)), ...
                this.t1w_on_pet);
            
            % clean
            if ~this.debug
                t4rb.deleteFourdfp(t4rb.theImages);
                t4rb.deleteFourdfp(t4rb.theImagesOp(:,1));
            end

            try
                if max(this.resolve_err(this.t1w_on_pet)) > 8
                    this = this.flirt_t1w_to_pet(no_clobber=opts.noclobber);
                end
            catch
            end

            popd(pwd0);
        end

        function this = resolve_delays(this)
        end
    end

    methods (Static)  % helpers
        function m = image_mass(obj)
            if ~isa(obj, 'mlfourd.ImagingContext2')
                obj = mlfourd.ImagingContext2(obj);
            end
            dV = voxelVolume(obj);
            ic1 = obj.thresh(0);
            sumDensities = dipsum(ic1);            
            m = sumDensities*dV;
        end

        function fn = json(obj)
            if ~isa(obj, 'mlfourd.ImagingContext2')
                obj = mlfourd.ImagingContext2(obj);
            end
            fn = obj.fqfp + ".json";
        end

        function jsonrecode(in, field, out)
            try
                str = struct(stackstr(3), field);
                jsonrecode(in, str, 'filenameNew', out);
            catch ME
                handwarning(ME)
                dispdbg(str)
                str = struct(stackstr(3), 'unknown field value');
                jsonrecode(in, str, 'filenameNew', out);
            end
        end

        function fn = mat(obj)
            if ~isa(obj, 'mlfourd.ImagingContext2')
                obj = mlfourd.ImagingContext2(obj);
            end
            fn = obj.fqfp + ".mat";
        end

        function fn = niigz(obj)
            if ~isa(obj, 'mlfourd.ImagingContext2')
                obj = mlfourd.ImagingContext2(obj);
            end
            fn = obj.fqfp + ".nii.gz";
        end

        function e = resolve_err(ic)
            j = jsondecode(fileread(strcat(ic.fqfp, '.json')));
            rerr = str2double(j.mlfourdfp_SimpleT4ResolveBuilder.cost_final.pairs_rotation_error.err); 
            terr = str2double(j.mlfourdfp_SimpleT4ResolveBuilder.cost_final.pairs_translation_error.err);
            e = [rerr, terr];
        end
    end

    %% PRIVATE

    properties
        atl_
        blur_
        pet_
        t1w_
        t1w_brain_
        t1w_dlicv_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
