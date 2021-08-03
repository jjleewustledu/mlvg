classdef VisionEval
    
    properties
        atlas
        feet1st_dyn
        head1st_dyn
        roi
        
        fdg
        co1
        co2
        ho1
        ho2
        oo1
        oo2        
    end
    
    methods (Static)
        function unpackCNDA(varargin)
            ip = inputParser;
            addRequired(ip, 'folder', @isfolder)
            addParameter(ip, 'f', 'sub-%n_%t_%d_%s', @ischar) % 20210219160921_FDG_Dynamic_8__
                                                       % filename (%a=antenna  (coil) number, 
                                                       %           %b=basename, 
                                                       %           %c=comments, 
                                                       %           %d=description, 
                                                       %           %e=echo number, 
                                                       %           %f=folder name, 
                                                       %           %i=ID of patient, 
                                                       %           %j=seriesInstanceUID, 
                                                       %           %k=studyInstanceUID, 
                                                       %           %m=manufacturer, 
                                                       %           %n=name of patient, 
                                                       %           %p=protocol, 
                                                       %           %r=instance number, 
                                                       %           %s=series number, 
                                                       %           %t=time, 
                                                       %           %u=acquisition number, 
                                                       %           %v=vendor, 
                                                       %           %x=study ID; 
                                                       %           %z=sequence name; default 'twilite')
            addParameter(ip, 'i', 'n') % ignore derived, localizer and 2D images (y/n, default n)
            addParameter(ip, 'o', pwd) % output directory (omit to save to input folder)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            [~,w] = mlbash('which dcm2niix');
            assert(~isempty(w))            
            [~,w] = mlbash('which pigz');
            if ~isempty(w)
                z = 'y';
            else
                z = 'n';
            end
            if ~isfolder(ipr.o)
                mkdir(ipr.o)
            end
            
            mlbash(sprintf('dcm2niix -f %s -i %s -o %s -z %s %s', ipr.f, ipr.i, ipr.o, z, ipr.folder))
        end
    end
    
    methods
        
        %%
        
        function call(this)
            this = this.resolveAll();
            
            first = {this.co1 this.oo1 this.ho1 this.fdg};
            second = {this.co2 this.oo2 this.ho2 this.oo1};
            tracer1 = {'CO' 'OO' 'HO' 'FDG'};
            tracer2 = {'CO' 'OO' 'HO' 'OO'};
            for n = 1:4
                first_tac = this.estimateTac(first{n});
                second_tac = this.estimateTac(second{n});
                this.plotComparison(first_tac.fourdfp.img, second_tac.fourdfp.img, tracer1{n}, tracer2{n})
            end
        end
        function ic_tac = estimateTac(this, ic)
            %ic_ta = ic.timeAveraged;
            %ic_ta.save();
            %[resolved,t4] = this.resolveToAtlas__(ic_ta);
            %resolved.fsleyes(this.atlas.fqfilename)
            %ic_on_atl = this.applyT4__(t4, ic);
            %ic_tac = ic_on_atl.volumeAveraged(this.roi);
            ic_tac = ic.volumeAveraged(this.roi);
        end
        function plotComparison(this, tac1, tac2, tracer1, varargin)
            if isempty(varargin)
                tracer2 = tracer1;
            else
                tracer2 = varargin{1};
            end
            
            taus1 = mlvg.SessionData.consoleTaus(tracer1);
            timesMid1 = cumsum(taus1) - taus1/2;
            taus2 = mlvg.SessionData.consoleTaus(tracer2);
            timesMid2 = cumsum(taus2) - taus2/2;
            
            h = figure; 
            plot(timesMid1, tac1, ':>', timesMid2, tac2, ':<')
            title([tracer1 ' activity in whole-brain'])
            legend('feet-1st', 'head-1st')
            xlabel('mid-frame time (s)')
            ylabel('activity (Bq/mL)')
            
            try
                savefig(h, ...
                    fullfile( ...
                    sprintf('VisionEval_plotComparison_%s_in_%s.fig', tracer1, this.roi.fileprefix)))
                figs = get(0, 'children');
                saveas(figs(1), ...
                    fullfile( ...
                    sprintf('VisionEval_plotComparison_%s_in_%s.png', tracer1, this.roi.fileprefix)))
                close(figs(1))
            catch ME
                handwarning(ME)
            end
        end
        function this = resolveAll(this)
            import mlfourd.ImagingContext2
            fdgs = globT('*_FDG_Dynamic_*.nii.gz');
            fdgs = fdgs(~contains(fdgs, '_avgt'));
            this.fdg = ImagingContext2(fdgs{1});
            this.fdg = this.fdg.selectFourdfp();
            this.fdg.fileprefix = 'fdg';
            fdg_avgt = copy(this.fdg);
            fdg_avgt = fdg_avgt.timeAveraged();
            fdg_early = copy(this.fdg.fourdfp);
            fdg_early.img = fdg_early.img(:,:,:,1:12);
            fdg_early = ImagingContext2(fdg_early);
            fdg_early_avgt = fdg_early.timeAveraged();
            fdg_early_avgt.fileprefix = 'fdg_early_avgt';
            
            cos = globT('*_CO_Dynamic_*.nii.gz');
            cos = cos(~contains(cos, '_avgt'));
            this.co1 = ImagingContext2(cos{1});
            this.co1 = this.co1.selectFourdfp();
            this.co1.fileprefix = 'co1';
            co1_avgt = this.co1.timeAveraged();
            this.co2 = ImagingContext2(cos{2});
            this.co2 = this.co2.selectFourdfp();
            this.co2.fileprefix = 'co2';
            co2_avgt = this.co2.timeAveraged();
            
            oos = globT('*_Oxygen_Dynamic_*.nii.gz');            
            oos = oos(~contains(oos, '_avgt'));
            this.oo1 = ImagingContext2(oos{1});
            this.oo1 = this.oo1.selectFourdfp();
            this.oo1.fileprefix = 'oo1';
            oo1_avgt = this.oo1.timeAveraged();
            this.oo2 = ImagingContext2(oos{2});
            this.oo2 = this.oo2.selectFourdfp();
            this.oo2.fileprefix = 'oo2';
            oo2_avgt = this.oo2.timeAveraged();
            
            hos = globT('*_Water_Dynamic_*.nii.gz');
            hos = hos(~contains(hos, '_avgt'));
            this.ho1 = ImagingContext2(hos{1});
            this.ho1 = this.ho1.selectFourdfp();
            this.ho1.fileprefix = 'ho1';
            ho1_avgt = this.ho1.timeAveraged();
            this.ho2 = ImagingContext2(hos{2});
            this.ho2 = this.ho2.selectFourdfp();
            this.ho2.fileprefix = 'ho2';
            ho2_avgt = this.ho2.timeAveraged();
            
            [~,t4_co_co] = this.resolveToTracer__(co2_avgt, co1_avgt);
            [~,t4_ho_ho] = this.resolveToTracer__(ho2_avgt, ho1_avgt);
            [~,t4_oo_oo] = this.resolveToTracer__(oo2_avgt, oo1_avgt);
            [~,t4_co_fdg] = this.resolveToTracer__(co1_avgt, fdg_early_avgt);
            [~,t4_ho_fdg] = this.resolveToTracer__(ho1_avgt, fdg_avgt);
            [~,t4_oo_ho] = this.resolveToTracer__(oo1_avgt, ho1_avgt);
            
            this.co1 = this.applyT4Composition__({         t4_co_fdg}, this.co1);
            this.co2 = this.applyT4Composition__({t4_co_co t4_co_fdg}, this.co2);
            this.ho1 = this.applyT4Composition__({         t4_ho_fdg}, this.ho1);
            this.ho2 = this.applyT4Composition__({t4_ho_ho t4_ho_fdg}, this.ho2);
            this.oo1 = this.applyT4Composition__({t4_oo_ho t4_ho_fdg}, this.oo1);
            this.oo2 = this.applyT4Composition__({t4_oo_oo t4_oo_ho t4_ho_fdg}, this.oo2);
        end
        function saveAll(this)
            this.co1.save
            this.co2.save
            this.ho1.save
            this.ho2.save
            this.oo1.save
            this.oo2.save
        end
        
        function this = VisionEval()           
            this.atlas = mlfourd.ImagingContext2( ...
                fullfile(getenv('REFDIR'), '711-2B_111.4dfp.hdr'));   
            this.atlas.filepath = pwd;
            this.atlas.save()
            this.roi = mlfourd.ImagingContext2( ...
                fullfile(getenv('REFDIR'), '711-2B_111_brain.4dfp.hdr')); 
            this.roi.filepath = pwd;
            this.roi.save()
        end
    end
    
    %% PRIVATE
    
    methods (Access = private)
        function out = applyT4__(~, t4, in)
            v = mlfourdfp.FourdfpVisitor;
            pwd0 = pushd(in.filepath);
            fp111 = '711-2B_111_op_111';
            if ~isfile(in.fqfilename)
                in.save()
            end
            out = v.t4img_4dfp(t4, in.fileprefix, ...
                               'out', [in.fileprefix '_op_111'], ...
                               'options', ['-O' fp111]);
            out = mlfourd.ImagingContext2([out '.4dfp.hdr']);
            popd(pwd0)
        end
        function out = applyT4Composition__(~, t4s, in)
            v = mlfourdfp.FourdfpVisitor;
            pwd0 = pushd(in.filepath);
            switch length(t4s)
                case 3                    
                    t4_ = v.t4_mul(t4s{1}, t4s{2});
                    t4 = v.t4_mul(t4_, t4s{3});
                case 2
                    t4 = v.t4_mul(t4s{1}, t4s{2});
                case 1
                    t4 = t4s{1};
                otherwise
                    error('mlvg:ValueError', ...
                        'VisionEval.applyT4Compostion:  length(t4) = %i', length(t4s))
            end
            fpfdg = 'fdg_avgt';
            if ~isfile(in.fqfilename)
                in.save()
            end
            out = v.t4img_4dfp(t4, in.fileprefix, ...
                               'options', ['-O' fpfdg]);
            out = mlfourd.ImagingContext2([out '.4dfp.hdr']);
            popd(pwd0)
        end
        function [resolved,t4] = resolveToAtlas__(this, ic3d)
            bldr = mlfourdfp.SimpleT4ResolveBuilder( ...
                'resolveTag', 'op_111', ...
                'theImages', {this.atlas, ic3d});
            bldr.resolve();
            resolved = mlfourd.ImagingContext2( ...
                [ic3d.fqfileprefix '_op_111.4dfp.hdr']);
            t4 = [ic3d.fqfileprefix '_to_op_111_t4'];
        end
        function [resolved,t4] = resolveToTracer__(~, source, target)
            bldr = mlfourdfp.SimpleT4ResolveBuilder( ...
                'resolveTag', ['op_' target.fileprefix], ...
                'theImages', {target, source});
            bldr.resolve();
            resolved = mlfourd.ImagingContext2( ...
                [source.fqfileprefix '_op_' target.fileprefix '.4dfp.hdr']);
            t4 = [source.fqfileprefix '_to_op_' target.fileprefix '_t4'];
        end
    end
end
