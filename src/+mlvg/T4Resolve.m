classdef T4Resolve < handle & mlsystem.IHandle
    %% run on machine.neuroimage.wustl.edu
    %  
    %  Created 04-Oct-2025 00:11:23 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.
    

    properties
        debug
    end

    properties (Dependent)
        unaligned
        unaligned_on_pet
        pet
        t1w
        t1w_brain
        t1w_dlicv
        t1w_on_pet
        t4rb
    end

    methods  % get, set
        function g = get.unaligned(this)
            g = this.unaligned_;
        end
        function g = get.unaligned_on_pet(this)
            assert(~isempty(this.unaligned))
            fqfn = fullfile(this.unaligned.filepath, strrep(this.unaligned.filename, 'unaligned', 'delay'));
            g = mlfourd.ImagingContext2(fqfn);
        end
        function g = get.pet(this)
            g = this.pet_;
            assert(isfile(g))  % overloaded by ImagingContext2
        end
        function g = get.t1w(this)
            assert(contains(this.t1w_.fileprefix, 'orient-std'))
            g = this.t1w_;
            assert(isfile(g))
        end
        function g = get.t1w_brain(this)
            if ~isempty(this.t1w_brain_)
                g = this.t1w_brain_;
                return
            end

            % read from filesystem
            fqfn = strcat(this.t1w.fqfp, '_brain.nii.gz');
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
            filepath = strrep(this.pet.filepath, 'sourcedata', 'derivatives');
            fqfn = fullfile(filepath, strcat('T1w_on_', this.pet.filename));
            g = mlfourd.ImagingContext2(fqfn);
        end
        function g = get.t4rb(this)
            g = this.t4rb_;
        end
    end

    methods
        function this = T4Resolve(opts)
            %% T4RESOLVE
            
            arguments
                opts.pet {mustBeNonempty}  % 3D representation, e.g. static, avgt, mipt
                opts.t1w = []
                opts.atl = fullfile(getenv("REFDIR"), "MNI152_T1_1mm.nii.gz")
                opts.blur double {mustBeScalarOrEmpty} = 6
                opts.debug logical = true
                opts.unaligned = []
            end
            
            this.pet_ = mlfourd.ImagingContext2(opts.pet);
            this.pet_.fqfn = convertStringsToChars(this.pet_.fqfn);
            this.t1w_ = mlfourd.ImagingContext2(opts.t1w);
            this.t1w_.fqfn = convertStringsToChars(this.t1w_.fqfn);
            assert(3 == ndims(this.pet_))
            this.atl_ = mlfourd.ImagingContext2(opts.atl);
            this.blur_ = opts.blur;
            this.debug = opts.debug;

            if ~isempty(opts.unaligned)
                this.unaligned_ = mlfourd.ImagingContext2(opts.unaligned);
            end

            %warning("off", "MATLAB:json:ExpectedNameOrEnd");
            %warning("off", "MATLAB:json:ExpectedCommaOrEnd");
        end

        function t1w_on_pet = find_t1w_on(this, pet)
            %% returns fqfp

            filepath = fileparts(pet);
            filepath = strrep(filepath, "sourcedata", "derivatives");
            t1w_on_pet = fullfile(filepath, ...
                sprintf('%s_op_%s', mybasename(this.t1w_brain), mybasename(pet)));
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
                'in', t1w__, ...
                'ref', this.pet, ...
                'out', this.t1w_on_pet, ...
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
        
        function this = resolve_unaligned_to_pet(this, opts)
            arguments
                this mlvg.T4Resolve
                opts.noclobber logical = false
            end

            if opts.noclobber && isfile(this.static_to_nmaf(this.unaligned_on_pet.fqfn))
                return
            end

            derivspath = strrep(this.pet.filepath, 'sourcedata', 'derivatives');
            pwd0 = pushd(derivspath);
            
            % resolve
            t1w_1 = mlfourd.ImagingContext2( ...
                fullfile(derivspath, strcat('T1w_on_', this.pet.filename)));
            msk_1 = t1w_1.binarized;
            msk_1.save();
            t1w_2 = mlfourd.ImagingContext2( ...
                fullfile(derivspath, strcat('T1w_on_', this.unaligned.filename)));
            msk_2 = t1w_2.binarized;
            msk_2.save();
            msks{1} = msk_1;
            msks{2} = msk_2;
            imgs{1} = this.pet;
            imgs{2} = this.unaligned;
            t4rb__ = mlfourd.SimpleT4ResolveBuilder( ...
                'workpath', derivspath, ...
                'maskForImages', msks, ...
                'theImages', imgs, ...
                'debug', this.debug);
            dispdbg(t4rb__)
            t4rb__ = t4rb__.resolve();
            
            % move files, write json
            copyfile(this.niigz(t4rb__.theImagesFinal{2}), this.niigz(this.unaligned_on_pet));
            delete(this.niigz(t4rb__.theImagesFinal{2}));
            j = this.json(t4rb__.theImagesFinal{2});
            if isfile(j)
                copyfile(j, this.json(this.unaligned_on_pet));
                delete(j);
            else
                copyfile(this.json(this.pet), this.json(this.unaligned_on_pet));
            end

            % t4img: nmaf unaligned -> nmaf delay0
            t4 = this.find_t4(this.unaligned, this.pet);
            unaligned_nmaf = myfileprefix(this.fourdfp(this.static_to_nmaf(this.unaligned.fqfn)));
            out = this.static_to_nmaf(this.unaligned_on_pet.fqfp);
            ref = this.static_to_nmaf(this.unaligned.fqfp);
            out1 = this.t4img_alt( ...
                t4, ...
                unaligned_nmaf, ...
                out=out, ...
                ref=ref);
            
            % write json of nmaf
            copyfile(this.json(unaligned_nmaf), this.json(out1));
            this.jsonrecode( ...
                out1, ...
                struct('image_activity', this.image_mass(out1)), ...
                out1);
            
            % clean
            if ~this.debug
                t4rb__.deleteFourdfp(t4rb__.theImages);
                t4rb__.deleteFourdfp(t4rb__.theImagesOp(:,1));
            end

            % store
            this.t4rb_ = t4rb__;

            popd(pwd0);
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
            t4rb__ = mlfourd.SimpleT4ResolveBuilder( ...
                'workpath', filepath, ...
                'maskForImages', msks, ...
                'theImages', imgs, ...
                'debug', this.debug);
            dispdbg(t4rb__)
            t4rb__ = t4rb__.resolve();
            
            % move files, write json
            movefile(this.niigz(t4rb__.theImagesFinal{2}), this.niigz(this.t1w_on_pet));
            j = this.json(t4rb__.theImagesFinal{2});
            if isfile(j)
                movefile(j, this.json(this.t1w_on_pet));
            else
                copyfile(this.json(this.t1w_brain), this.json(this.t1w_on_pet));
            end
            this.jsonrecode( ...
                this.t1w_on_pet, ...
                struct('image_activity', this.image_mass(this.t1w_on_pet)), ...
                this.t1w_on_pet);
            
            % clean
            if ~this.debug
                t4rb__.deleteFourdfp(t4rb__.theImages);
                t4rb__.deleteFourdfp(t4rb__.theImagesOp(:,1));
            end

            % store
            this.t4rb_ = t4rb__;

            try
                if max(this.resolve_err(this.t1w_on_pet)) > 8
                    this = this.flirt_t1w_to_pet(no_clobber=opts.noclobber);
                end
            catch
            end

            popd(pwd0);
        end
        
        function this = resolve_t1w_to_intermed_to_pet(this, opts)
            arguments
                this mlvg.T4Resolve
                opts.intermed {mustBeNonempty}
                opts.noclobber logical = false
            end

            if opts.noclobber && isfile(this.t1w_on_pet)
                return
            end

            filepath = strrep(this.pet.filepath, 'sourcedata', 'derivatives');
            pwd0 = pushd(filepath);
            
            % resolve intermed to pet
            intermed = mlfourd.ImagingContext2(char(opts.intermed));
            intermed_dlicv = this.make_dlicv(intermed);
            msks{1} = mlfourd.ImagingContext2('none.nii.gz');
            msks{2} = intermed_dlicv;
            imgs{1} = this.pet;
            imgs{2} = intermed;
            t4rb__ = mlfourd.SimpleT4ResolveBuilder( ...
                'workpath', filepath, ...
                'maskForImages', msks, ...
                'theImages', imgs, ...
                'debug', this.debug);
            dispdbg(t4rb__)
            t4rb__ = t4rb__.resolve();

            % t4img: t1w -> pet
            t4 = this.find_t4(intermed, this.pet);
            t1w_on_intermed = this.find_t1w_on(intermed);
            ref = fullfile(filepath, this.pet.fileprefix);
            out = fullfile(filepath, this.t1w_on_pet.fileprefix);
            out1 = this.t4img( ...
                t4, ...
                t1w_on_intermed, ...
                out=out, ...
                ref=ref);
            
            % write json
            j = this.json(t4rb__.theImagesFinal{2});
            if isfile(j)
                movefile(j, this.json(out1));
            else
                copyfile(this.json(this.t1w_brain), this.json(out1));
            end
            this.jsonrecode( ...
                out1, ...
                struct('image_activity', this.image_mass(out1)), ...
                out1);
            
            % clean
            if ~this.debug
                t4rb__.deleteFourdfp(t4rb__.theImages);
                t4rb__.deleteFourdfp(t4rb__.theImagesOp(:,1));
            end

            % store
            this.t4rb_ = t4rb__;

            popd(pwd0);
        end

        function this = resolve_t1w_to_intermeds_to_pet(this, opts)
            arguments
                this mlvg.T4Resolve
                opts.intermed {mustBeNonempty}
                opts.intermed2 {mustBeNonempty}
                opts.noclobber logical = false
            end

            if opts.noclobber && isfile(this.t1w_on_pet)
                return
            end

            filepath = strrep(this.pet.filepath, 'sourcedata', 'derivatives');
            pwd0 = pushd(filepath);
            
            % resolve intermed to pet
            intermed = mlfourd.ImagingContext2(char(opts.intermed));
            % intermed_dlicv = this.make_dlicv(intermed);
            intermed2 = mlfourd.ImagingContext2(char(opts.intermed2));
            msks{1} = mlfourd.ImagingContext2('none.nii.gz');
            msks{2} = mlfourd.ImagingContext2('none.nii.gz');  % intermed_dlicv;
            imgs{1} = this.pet;
            imgs{2} = intermed2;
            t4rb__ = mlfourd.SimpleT4ResolveBuilder( ...
                'workpath', filepath, ...
                'maskForImages', msks, ...
                'theImages', imgs, ...
                'debug', this.debug);
            dispdbg(t4rb__)
            t4rb__ = t4rb__.resolve();

            % t4img: t1w -> pet
            t4 = this.find_t4(intermed2, this.pet);
            t1w_on_intermed = this.find_t1w_on(intermed);
            ref = fullfile(filepath, this.pet.fileprefix);
            out = fullfile(filepath, this.t1w_on_pet.fileprefix);
            out1 = this.t4img( ...
                t4, ...
                t1w_on_intermed, ...
                out=out, ...
                ref=ref);
            
            % write json
            j = this.json(t4rb__.theImagesFinal{2});
            if isfile(j)
                movefile(j, this.json(out1));
            else
                copyfile(this.json(this.t1w_brain), this.json(out1));
            end
            this.jsonrecode( ...
                out1, ...
                struct('image_activity', this.image_mass(out1)), ...
                out1);
            
            % clean
            if ~this.debug
                t4rb__.deleteFourdfp(t4rb__.theImages);
                t4rb__.deleteFourdfp(t4rb__.theImagesOp(:,1));
            end

            % store
            this.t4rb_ = t4rb__;

            % try simpler alternatives when resolve_err is high
            try
                if max(this.resolve_err(this.t1w_on_pet)) > 3.5
                    this = this.resolve_t1w_to_pet(noclobber=opts.noclobber);
                end
            catch
            end

            popd(pwd0);
        end

        function intermed_dlicv = make_dlicv(this, intermed)
            arguments
                this mlvg.T4Resolve
                intermed mlfourd.ImagingContext2
            end

            t1wb_on_intermed = fullfile( ...
                fileparts(intermed), ...
                strcat('T1w_on_', intermed.filename));
            assert(isfile(t1wb_on_intermed))
            t1wb_on_intermed = mlfourd.ImagingContext2(t1wb_on_intermed);
            intermed_dlicv = t1wb_on_intermed.binarized();
            intermed_dlicv.fileprefix = strcat(intermed.fileprefix, '_DLICV');
            intermed_dlicv.save();
        end

        function this = resolve_delays(this)
        end
    end

    methods (Static)  % helpers
        function t4 = find_t4(src, dest)
            arguments
                src {mustBeNonempty}
                dest {mustBeNonempty}
            end
            if isa(src, "mlio.IOInterface")
                src = src.fqfn;
            end
            if isa(dest, "mlio.IOInterface")
                dest = dest.fqfn;
            end

            filepath = fileparts(dest);
            filepath = strrep(filepath, 'sourcedata', 'derivatives');
            filepath = regexprep(filepath, 'ses-\d{14}', 'ses-*');
            bsrc = mybasename(src);
            bdest = mybasename(dest);
            t4 = mglob(fullfile(filepath, 'Log', sprintf('%s_to_op_%s_t4', bsrc, bdest)), aschar=true);
            if isemptytext(t4)
                t4 = mglob(fullfile(filepath, sprintf('%s_to_op_%s_t4', bsrc, bdest)), aschar=true);
            end
            if iscell(t4)
                t4 = t4{1};
            end
            assert(ischar(t4))  % not cell, not string
        end

        function fqfn = fourdfp(obj, opts)
            %% manages faulty center coords of 4dfp

            arguments
                obj {mustBeNonempty}
                opts.ref = []
            end
            ic = mlfourd.ImagingContext2(obj);
            if contains(ic.filesuffix, '.4dfp') % return trivial
                fqfn = ic.fqfn;
                return
            end

            ic.selectFourdfpTool();
            ifc = ic.fourdfp();
            if ~isempty(opts.ref)
                ifc_ref = mlfourd.ImagingFormatContext2(opts.ref);
                assert(contains(ifc_ref.filesuffix, '.4dfp'))
                ifc.hdr = ifc_ref.hdr;
            end
            ifc.save();
            fqfn = ifc.fqfn;
        end 

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

        function fqfn = niigz(obj, opts)
            %% manages faulty center coords of 4dfp; returns fqfn

            arguments
                obj {mustBeNonempty}
                opts.ref = []
            end
            ic = mlfourd.ImagingContext2(obj);
            if contains(ic.filesuffix, '.nii') % return trivial
                fqfn = ic.fqfn;
                return
            end

            ic.selectFourdfpTool();
            ifc = ic.nifti();
            if ~isempty(opts.ref)
                ifc_ref = mlfourd.ImagingFormatContext2(opts.ref);
                assert(contains(ifc_ref.filesuffix, '.nii'))
                ifc.hdr = ifc_ref.hdr;
            end
            ifc.save();
            fqfn = ifc.fqfn;
        end 

        function e = resolve_err(ic)
            j = jsondecode(fileread(strcat(ic.fqfp, '.json')));
            rerr = str2double(j.mlfourdfp_SimpleT4ResolveBuilder.cost_final.pairs_rotation_error.err); 
            terr = str2double(j.mlfourdfp_SimpleT4ResolveBuilder.cost_final.pairs_translation_error.err);
            e = [rerr, terr];
        end

        function fqfp = static_to_nmaf(obj)
            if isa(obj, "mlfourd.ImagingContext2") || isa(obj, "mlfourd.ImagingFormatContext2")
                obj = obj.fqfp;
            end
            fqfp = strrep(obj, "createNiftiStatic", "createNiftiMovingAvgFrames");
            fqfp = char(fqfp);
        end

        function out = t4img(t4, in, opts)
            %% returns nii.gz of transformation, with due diligence with hdr

            arguments
                t4 {mustBeFile}
                in {mustBeNonempty}
                opts.out {mustBeNonempty}
                opts.ref {mustBeNonempty}
            end

            pwd0 = pushd(myfileparts(in));
            bv = mlfourdfp.FourdfpVisitor();
            bv.t4img_4dfp( ...
                t4, ...
                strcat(myfileprefix(in), '.4dfp.img'), ...
                out=myfileprefix(opts.out), ...
                options=strcat('-O', myfileprefix(opts.ref)));
            out = mlvg.T4Resolve.niigz( ...
                strcat(myfileprefix(opts.out), '.4dfp.img'), ...
                ref=strcat(myfileprefix(opts.ref), '.nii.gz'));
            popd(pwd0);
        end

        function out = t4img_alt(t4, in, opts)
            %% directly make system calls to t4img_4dfp;
            % returns nii.gz of transformation, with due diligence with hdr

            arguments
                t4 {mustBeFile}
                in {mustBeTextScalar}
                opts.out {mustBeTextScalar}
                opts.ref {mustBeTextScalar}
            end
            in = myfileprefix(in);
            opts.out = myfileprefix(opts.out);
            opts.ref = myfileprefix(opts.ref);

            pwd0 = pushd(myfileparts(in));
            cmd = sprintf('t4img_4dfp %s %s %s -O%s', t4, in, opts.out, opts.ref);
            mysystem(cmd);
            out = mlvg.T4Resolve.niigz( ...
                strcat(myfileprefix(opts.out), '.4dfp.img'), ...
                ref=strcat(myfileprefix(opts.ref), '.nii.gz'));
            popd(pwd0);
        end
    end

    %% PRIVATE

    properties (Access = private)
        atl_
        blur_
        unaligned_
        pet_
        t1w_
        t1w_brain_
        t1w_dlicv_
        t4rb_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
