classdef TracerSuvrBuilder < mlfourdfp.AbstractSessionBuilder
	%% TRACERSUVRBUILDER works on a single subject at a time

	%  $Revision$
 	%  was created 28-Mar-2018 22:00:52 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlpet/src/+mlpet.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.  Copyright 2018 John Joowon Lee.

    properties (Constant)
        NORMAL_CBF = 44 / 1.05 % mL/hg/min, density of brain ~ 1.05 g/mL
        NORMAL_CBV = 3.8 / 1.05 % mL/hg
        NORMAL_CMRGLC = 30.0; % umol/hg/min
        NORMAL_CMRO2 = 3.3 / 1.05 % mL/hg/min
        NORMAL_OEF = 0.44
        SUPPORTED_TRACERS = {'HO' 'OO' 'OC'} % 1st is ReferenceTracer
    end
    
    properties 
        atlasVoxelSize = 222;
        outpath
        rebuild = false;
        tracerKind = 'tracerResamplingRestricted' % method@SessionData
        workpath
    end
    
	properties (Dependent)  
        ReferenceTracer
        supportedTracers
    end
    
    methods (Static)        
        function prod = averageProduct(prod)
            assert(iscell(prod))
            p1 = prod{1};
            assert(isa(p1, 'mlpet.SuvrContext'))
            for idx = 2:length(prod)
                p1 = p1 + prod{idx};
            end
            prod = p1 ./ length(prod);
            pos = regexp(prod.fileprefix, '_on_T1001');
            pos = pos - 1;
            prod.fileprefix(pos-5:pos) = '000000';
        end
        function globbed = globTracer(tr)
            assert(ischar(tr))
            tr = lower(tr);
            globbed_ = globT([tr 'dt*.4dfp.hdr']);
            globbed = {};
            for g = globbed_
                if ~lstrfind(g{1}, '_avgt')
                    globbed = [globbed g{1}]; %#ok<*AGROW>
                end
            end
        end
        function suvrCon = t4img_to_T1(suvrCon)
            %% works in pwd
            
            fv = mlfourdfp.FourdfpVisitor();
            suvrCon.save
            ss = strsplit(suvrCon.fileprefix, '_times');
            t4 = sprintf('%s_to_T1001_t4', ss{1});
            targ = sprintf('%s_on_T1001', suvrCon.fileprefix);
            fv.t4img_4dfp(t4, suvrCon.fileprefix, 'out', targ, 'options', '-OT1001')
            suvrCon = mlpet.SuvrContext('sessionData', suvrCon.sessionData, ...
                                        'filename', [targ '.4dfp.hdr']);
        end
    end

	methods   
        
        %% GET/SET
        
        function g = get.ReferenceTracer(this)
            g = this.SUPPORTED_TRACERS{1};
        end
        function g = get.supportedTracers(this)
            g = lower(this.SUPPORTED_TRACERS);
        end
        
        %%
        
        function this = buildAll(this)
            %% top-level build
            %  @returns this.product := {cbf cbv y cmro2 oef mask} in physiol. units
            
            pwd0 = pushd(this.workpath);
            
            theCbf = this.averageProduct(this.buildCbf.product);
            theCbv = this.averageProduct(this.buildCbv.product);
            theY   = this.averageProduct(this.buildY.product);
            this   = this.buildBetas(theCbf, theCbv, theY);
            this.product_ = [{theCbf  theCbv theY} this.product];      
            
            popd(pwd0)
        end
        
        
        
        
        function this = buildCbf(this)
            this = this.buildTracer('tracer', 'ho', 'physiol', 'cbf', 'expected', this.NORMAL_CBF);
        end
        function this = buildCbv(this)
            this = this.buildTracer('tracer', 'oc', 'physiol', 'cbv', 'expected', this.NORMAL_CBV);
        end
        function this = buildCmrglc(this)
            this = this.buildTracer('tracer', 'fdg', 'physiol', 'cmrglc', 'expected', this.NORMAL_CMRGLC);
        end
        function [this,ogi] = buildGlcMetab(this)
            assert(lexist(this.tracerSuvrNamed('cmro2'), 'file'));            
            assert(lexist(this.tracerSuvrNamed('fdg'),   'file'));
            cmro2  = this.tracerSuvrNamed('cmro2', 'typ', 'numericalNiftid');
            cmrglc = this.tracerSuvrNamed('fdg',   'typ', 'numericalNiftid');
            
            ogi = (cmro2 ./ cmrglc) * 5.4;
            ogi.fqfilename = this.tracerSuvrNamed('ogi');
            ogi.save;
            ogi = mlfourd.ImagingContext2(ogi);
            
            this.product_ = ogi;
        end 
        function this = buildY(this)
            this = this.buildTracer('tracer', 'oo', 'physiol', 'Y');
        end
        function this = buildTracer(this, varargin)
            ip = inputParser;
            addParameter(ip, 'tracer', 'ho')
            addParameter(ip, 'physiol', 'cbf')
            addParameter(ip, 'expected', 1, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            pwd0 = pushd(this.workpath);
            
            prods = {};
            tracs = this.globTracer(ipr.tracer);
            for t = tracs
                sc = mlpet.SuvrContext('sessionData', this.sessionData, 'filename', t{1});
                sc = sc.timeAveraged('suffix', 'none');
                sc.fqfilename = sc.fqfilenameTimeWindowed();
                sc = this.t4img_to_T1(sc);
                sc = this.convertToPhysiol(sc, ipr.expected);
                sc = sc.blurred(4.3);
                sc.fqfileprefix = fullfile(this.outpath, [ipr.physiol 'dt' sc.datestr() '_on_T1001']);
                sc.save
                this.saveasNifti(sc)
                mlbash(sprintf('mv -f *_times*.4dfp.* %s', this.outpath))
                prods = [prods {sc}];
            end
            this.product_ = prods;
            
            popd(pwd0)
        end        
        function [this,mdl] = buildBetas(this, cbf, cbv, y)
            msk  = this.constructMaskContext();
            msk_ = logical(msk.fourdfp.img);
            cbf_ = ensureColVector(squeeze(cbf.fourdfp.img(msk_))) / this.NORMAL_CBF;
            cbv_ = ensureColVector(squeeze(cbv.fourdfp.img(msk_))) / this.NORMAL_CBV;
            y_   = ensureColVector(squeeze(  y.fourdfp.img(msk_)));
            
            % nonlinear regression
            tbl = table(cbf_, cbv_, y_);            
            mdlfun = @(b,x) b(1)*x(:,1) + b(2)*x(:,2);
            beta0 = [1 1];
            mdl = fitnlm(tbl, mdlfun, beta0);
            beta1 = mdl.Coefficients{1, 'Estimate'};
            beta2 = mdl.Coefficients{2, 'Estimate'};
            disp(mdl)           
            fprintf('mdl.RMSE->%g \t min(y_)->%g \t max(y_)->%g\n', mdl.RMSE, min(y_), max(y_))
            
            % assign cmro2, oef
            cmro2 = y - cbv .* beta2;
            cmro2 = cmro2.scrubNegative();
            cmro2 = this.convertToPhysiol(cmro2, this.NORMAL_CMRO2);
            cmro2.fqfilename = fullfile(this.outpath, ['cmro2dt' this.datestr() '_on_T1001']);
            cmro2.save;  
            this.saveasNifti(cmro2)
            
            oef = cmro2 ./ (cbf .* beta1);
            oef = oef .* msk;
            oef = oef.scrubNanInf();
            oef = oef.scrubNegative();
            oef = this.convertToPhysiol(oef, this.NORMAL_OEF);
            oef = oef.scrubGt(1);
            oef.fqfilename = fullfile(this.outpath, ['oefdt' this.datestr() '_on_T1001']);
            oef.save; 
            this.saveasNifti(oef)
                   
            this.product_ = {cmro2 oef msk};
        end  
                
        %% Utilities        
        
        function ic = constructMaskContext(this)
            ic = mlfourd.ImagingContext2(fullfile(this.workpath, 'wmparc.4dfp.hdr'));            
            ic = ic.blurred(4.3);
            ic = ic.binarized;
        end
        function con = convertToPhysiol(this, con, expected)
            assert(isa(con, 'mlpet.SuvrContext'))
            assert(isscalar(expected))
            volavg = con.maskedMean(this.constructMaskContext());
            con = con .* expected ./ volavg;
        end
        function d = datestr(this)
            d = datestr(this.sessionData.datetime, 'yyyymmddHHMMSS');
            d = [d(1:8) '000000'];
        end
        
        %%
        
 		function this = TracerSuvrBuilder(varargin)
 			this = this@mlfourdfp.AbstractSessionBuilder(varargin{:});
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'workpath', this.sessionData.dataPath)
            addParameter(ip, 'outpath',  fullfile(this.sessionData.dataPath, 'SUVR', ''))
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.workpath = ipr.workpath;
            this.outpath = ipr.outpath;
            ensuredir(this.outpath)
 		end
    end 
    
    %% PRIVATE
    
    properties (Access = private)
    end
    
    methods (Access = private)
    end
    
    %% HIDDEN
    
    methods (Hidden)  
        function saveasNifti(~, suvrCon)
            suvrCon.filesuffix = '.nii.gz';
            suvrCon.save
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

