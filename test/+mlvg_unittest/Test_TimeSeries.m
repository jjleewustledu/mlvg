classdef Test_TimeSeries < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 10-Jul-2024 15:41:01 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/test/+mlvg_unittest.
    %  Developed on Matlab 24.1.0.2628055 (R2024a) Update 4 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        pet = "sub-108237_ses-20221031100910_trc-co_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames-ParcSchaeffer-reshape-to-schaeffer-schaeffer.nii.gz"
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mlvg.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
        function test_analytic_signal_arousal_coarse(this)    

            statistic = "median";
            % plot_kind = "abs";
            % plot_kind = "angle";
            plot_kind = "unwrap_angle";
            hp_thresh_hz = 0.01;
            lp_thresh_hz = 0.05;
            decimate = 3;
            if contains(plot_kind, "angle")
                cmap = @viridis;
            else
                cmap = @magma;
            end

            pwd0 = pushd("/Volumes/PrecunealSSD/Singularity/CCIR_01211/derivatives");
            for t = 3 % [10,5,3]
                figure
                tiledlayout(2, 3)
                globbed = asrow(mglob("sub-*/ses-*/pet/sub-*_ses-*_trc-co_proc-delay0-BrainMoCo2*tau"+t+"*reshape-to-schaeffer-schaeffer.nii.gz"));
                for g = flip(globbed)

                    obj = mlvg.TimeSeries(g, t0=60, decimate=decimate, hp_thresh_hz=hp_thresh_hz, lp_thresh_hz=lp_thresh_hz);
                    Nt = size(obj, 2);

                    % rois = ["skull", "arousal", "extrinsic", "dmn"];
                    rois = ["extrinsic", "dmn"];
                    % rois = ["cerebellar_cortex", "cerebellar_white", "vermis"];
                    % rois = ["vis", "sms", "dan", "van", "limbic", "con", "dmn"];
                    % rois = ["hippocampus", "amygdala", "striatum", "thalamus", "dmn"];
                    parc_img = nan(length(rois), Nt);
                    for ridx = 1:length(rois)
                        parc_img(ridx,:) = obj.parc(rois(ridx), statistic=statistic);
                    end

                    parc_img_nobolus = obj.build_centered_and_rescaled(parc_img');  % times x parcs
                    em = hilbert(parc_img_nobolus);

                    ref_img_nobolus = obj.build_centered_and_rescaled(obj.parc("arousal", statistic=statistic)');  % times x parcs
                    arousal = hilbert(ref_img_nobolus);

                    % as = obj.build_band_passed( ...
                    %     obj.build_centered_and_rescaled(conj(arousal) .* em));  % times x parcs
                    as = obj.build_band_passed(em) .* obj.build_band_passed(conj(arousal));  % times x parcs

                    nexttile;
                    times = ascol(0:decimate:decimate*Nt-1);
                    switch plot_kind
                        case "abs"
                            plot(times, abs(as), LineWidth=1.5);
                        case "angle"
                            plot(times, angle(as), LineWidth=1.5);
                            ylim([0, pi]);
                        case "unwrap_angle"
                            plot(times, unwrap(angle(as)), LineWidth=1.5);
                        otherwise
                            error("mlvg:ValueError", "%s: plot_kind->%s", stackstr(), plot_kind)
                    end
                    colororder(cmap(size(parc_img, 1)));
                    % colormap(cmap); colorbar;
                    legend(rois, Interpreter="none");
                    title(obj.fileprefix, FontSize=5, Interpreter="none")
                    if strcmp(plot_kind, "abs")
                        ylabel("amplitude (normalized)", FontSize=14)
                    else
                        ylabel("phase (radians)", FontSize=14)
                    end
                    xlabel("time excluding bolus passage (sec)", FontSize=14)

                    % xlim([0, 100]);
                    % ylim([-pi pi]);
                end
            end

            popd(pwd0);
        end
        function test_analytic_signal_arousal_all(this)

            % plot_kind = "abs";
            plot_kind = "unwrap_angle";
            % band_passed = true;
            band_passed = false;

            obj = this.testObj;
            cd("/Volumes/PrecunealSSD/Singularity/CCIR_01211/derivatives")
            for t = [10,5,3]
                figure
                tiledlayout(1, 6)
                globbed = asrow(mglob("sub-*/ses-*/pet/sub-*_ses-*_trc-co_proc-delay0-BrainMoCo2*tau"+t+"*reshape-to-schaeffer-schaeffer.nii.gz"));
                for g = globbed

                    co = mlfourd.ImagingFormatContext2(g);

                    co_img_nobolus = co.img(:, 61:end)';  % times x parcs
                    co_img_nobolus = obj.build_global_signal_regressed(co_img_nobolus);
                    em  = hilbert(co_img_nobolus);
                    arousal = hilbert(median(co_img_nobolus(:, [11, 17, 32]), 2));  % 4th ventricle, ventral dicencephalon
                    if band_passed
                        as = obj.build_band_passed( ...
                            obj.build_centered_and_rescaled(conj(arousal)) .* ...
                            obj.build_centered_and_rescaled(em));  % times x parcs
                    else
                        as = obj.build_centered_and_rescaled(conj(arousal)) .* ...
                            obj.build_centered_and_rescaled(em);  % times x parcs
                    end

                    nexttile;
                    switch plot_kind
                        case "abs"
                            plot(abs(as));
                        case "unwrap_angle"
                            plot(unwrap(angle(as)));
                        otherwise
                            error("mlvg:ValueError", "%s: plot_kind->%s", stackstr(), plot_kind)
                    end
                    colororder(viridis);
                    colormap(viridis); colorbar;
                    title(co.fileprefix, FontSize=6, Interpreter="none")
                end
            end
        end
        function test_analytic_signal_arousal(this)
            obj = this.testObj;
            cd("/Users/jjlee/Documents/PapersMine/PapersInProgress/PET_connectivity/CO")
            co = mlfourd.ImagingFormatContext2("sub-108237_ses-20221031100910_trc-co_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames-ParcSchaeffer-reshape-to-schaeffer-schaeffer.nii.gz");
            
            co_img_nobolus = co.img(:, 61:end)';  % times x parcs
            co_img_nobolus = obj.build_global_signal_regressed(co_img_nobolus);
            em  = hilbert(co_img_nobolus);
            arousal = hilbert(median(co_img_nobolus(:, [11, 17, 32]), 2));  % 4th ventricle, ventral dicencephalon
            % obj.build_band_passed( ...
            as = obj.build_centered_and_rescaled(conj(arousal)) .* ...
                obj.build_centered_and_rescaled(em);  % times x parcs
             
            figure; 
            plot(abs(as));
            colororder(viridis); 
            colormap(viridis); colorbar;
            
            figure; 
            plot(unwrap(angle(as)));
            colororder(viridis); 
            colormap(viridis); colorbar;

            % figure; plot(median(abs(as), 1));
            % figure; plot(median(unwrap(angle(as)), 1));

        end
        function test_analytic_signal_arbitrary(this)
            obj = this.testObj;
            cd("/Users/jjlee/Documents/PapersMine/PapersInProgress/PET_connectivity/CO")
            co = mlfourd.ImagingFormatContext2("sub-108237_ses-20221031100910_trc-co_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames-ParcSchaeffer-reshape-to-schaeffer-schaeffer.nii.gz");
            
            co_img_nobolus = co.img(:, 61:end)';  % times x parcs
            em  = hilbert(obj.build_global_signal_regressed(co_img_nobolus));
            arousal_indices = 41;  % head extracerebral
            %  130;  % air
            %  42:109;  % wm
            % [11, 17, 32]  % 4th ventricle, ventral dicencephalon
            arousal = hilbert(obj.build_centered_and_rescaled( ...
                median(co_img_nobolus(:, arousal_indices), 2)));
            as = obj.build_band_passed( ...
                obj.build_centered_and_rescaled(conj(arousal)) .* ...
                obj.build_centered_and_rescaled(em));  % times x parcs
             
            figure; 
            plot(abs(as));
            colororder(viridis); 
            colormap(viridis); colorbar;

            figure; 
            plot(unwrap(angle(as)));
            colororder(viridis); 
            colormap(viridis); colorbar;

            % figure; plot(median(abs(as), 1));
            % figure; plot(median(unwrap(angle(as)), 1));

        end
        function test_corrcoef_ori(this)
            cd("/Users/jjlee/Documents/PapersMine/PapersInProgress/PET_connectivity/CO")

            co = mlfourd.ImagingFormatContext2("sub-108237_ses-20221031100910_trc-co_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames-ParcSchaeffer-reshape-to-schaeffer-schaeffer.nii.gz");
            
            bolus = mean(co.img, 1);
            co_img_nobolus = co.img - asrow(bolus);  % remove mean bolus
            figure; imagesc(co_img_nobolus)
            cc = corrcoef(co_img_nobolus');
            figure; imagesc(cc)
        end
        function test_corrcoef_all(this)
            obj = this.testObj;
            cd("/Volumes/PrecunealSSD/Singularity/CCIR_01211/derivatives")
            for t = [10,5,3]
                figure
                tiledlayout(1, 6)
                globbed = asrow(mglob("sub-*/ses-*/pet/sub-*_ses-*_trc-co_proc-delay0-BrainMoCo2*tau"+t+"*reshape-to-schaeffer-schaeffer.nii.gz"));
                for g = globbed
                    co = mlfourd.ImagingFormatContext2(g);

                    co_img_nobolus = co.img(:, 61:end)';  % times x parcs
                    co_img_nobolus = ...
                        obj.build_band_passed( ...
                        obj.build_centered_and_rescaled(co_img_nobolus));
                    % co_img_nobolus = ...
                    %     obj.build_band_passed( ...
                    %     obj.build_centered_and_rescaled( ...
                    %     obj.build_global_signal_regressed(co_img_nobolus)));

                    nexttile; 
                    cc = corrcoef(co_img_nobolus);
                    imagesc(cc);
                    title(co.fileprefix, FontSize=6, Interpreter="none")
                end
            end
        end
        function test_corrcoef_single(this)
            obj = this.testObj;
            cd("/Users/jjlee/Documents/PapersMine/PapersInProgress/PET_connectivity/CO")
            co = mlfourd.ImagingFormatContext2("sub-108237_ses-20221031100910_trc-co_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames-ParcSchaeffer-reshape-to-schaeffer-schaeffer.nii.gz");
            
            % co_img_nobolus = co.img(:, 61:end)';  % times x parcs
            % co_img_nobolus = ...
            %     obj.build_band_passed( ...
            %     obj.build_centered_and_rescaled( ...
            %     obj.build_global_signal_regressed(co_img_nobolus)));

            co_img_nobolus = co.img;  % times x parcs
            co_img_nobolus = ...
                obj.build_band_passed( ...
                obj.build_centered_and_rescaled(co_img_nobolus));

            figure; imagesc(co_img_nobolus)
            cc = corrcoef(co_img_nobolus);
            figure; imagesc(cc)
        end
        function test_corrcoef(this)
            %% build correlation matrices for all CO schaeffer dynamic data, all frame window lengths

            for tau = [10,5,3]
                globbed = asrow(mglob(sprintf( ...
                    '/Users/jjlee/Singularity/CCIR_01211/derivatives/sub-*/ses-*/pet/sub-*_ses-*_trc-co_proc-*-tau%i-*-reshape-to-schaeffer-schaeffer.nii.gz', tau)));
                for g = globbed
                    try
                        pwd0 = pushd(fileparts(g));  % sub-*/ses-*/

                        co = mlfourd.ImagingFormatContext2(g);
                        % bolus = mean(co.img, 1);
                        % co_img_nobolus = co.img - asrow(bolus);  % remove mean bolus
                        co_img_nobolus = co.img(:,61:end);  % drop boluses in first 60 sec
                        % figure; 
                        % imagesc(co_img_nobolus)
                        % cc = corrcoef(co_img_nobolus(:, 35:end)');
                        cc = corrcoef(co_img_nobolus');
                        h = figure; 
                        imagesc(cc); colorbar;

                        fig_fp = myfileprefix(strrep(g, "-ParcSchaeffer-reshape-to-schaeffer-schaeffer", "-corrcoef"));
                        title(mybasename(fig_fp), Interpreter="none");
                        saveFigure2(h, fig_fp);

                        popd(pwd0);
                    catch ME
                        handwarning(ME)
                    end
                end
            end
        end
        function test_homotopy(this)
            [w,cw] = this.testObj.parc("precuneus_white");
            [g,cg] = this.testObj.parc("precuneus_grey");

            tiledlayout(2, 1) %, TileSpacing="compact", Padding="compact")
            nexttile
            this.testObj.plot(w, id="precuneus white", combo=cw, use_cbrewer2=false, colormap="cividis", LineWidth=2)
            nexttile
            this.testObj.plot(g, id="precuneus grey", combo=cg, use_cbrewer2=false, colormap="viridis", LineWidth=2)
            fontsize(scale=2);
        end
        function test_7rsns(this)
            [r1,c1] = this.testObj.parc("vis");
            [r2,c2] = this.testObj.parc("sms");
            [r3,c3] = this.testObj.parc("dan");
            [r4,c4] = this.testObj.parc("van");
            [r5,c5] = this.testObj.parc("limbic");
            [r6,c6] = this.testObj.parc("con");
            [r7,c7] = this.testObj.parc("dmn");

            tiledlayout(7, 1) %, TileSpacing="compact", Padding="compact")
            nexttile
            this.testObj.plot(r1, id="vis", combo=c1, use_cbrewer2=false, colormap="viridis")
            nexttile
            this.testObj.plot(r2, id="sms", combo=c2, use_cbrewer2=false, colormap="viridis")
            nexttile
            this.testObj.plot(r3, id="dan", combo=c3, use_cbrewer2=false, colormap="viridis")
            nexttile
            this.testObj.plot(r4, id="van", combo=c4, use_cbrewer2=false, colormap="viridis")
            nexttile
            this.testObj.plot(r5, id="limbic", combo=c5, use_cbrewer2=false, colormap="viridis")
            nexttile
            this.testObj.plot(r6, id="con", combo=c6, use_cbrewer2=false, colormap="viridis")
            nexttile
            this.testObj.plot(r7, id="dmn", combo=c7, use_cbrewer2=false, colormap="viridis")
            fontsize(scale=1.25);
        end
        function test_7rsns_mean(this)
            [r1,c1] = this.testObj.parc("vis", statistic="mean");
            [r2,c2] = this.testObj.parc("sms", statistic="mean");
            [r3,c3] = this.testObj.parc("dan", statistic="mean");
            [r4,c4] = this.testObj.parc("van", statistic="mean");
            [r5,c5] = this.testObj.parc("limbic", statistic="mean");
            [r6,c6] = this.testObj.parc("con", statistic="mean");
            [r7,c7] = this.testObj.parc("dmn", statistic="mean");

            tiledlayout(7, 1) %, TileSpacing="compact", Padding="compact")
            nexttile
            this.testObj.plot(r1, id="vis", combo=c1, use_cbrewer2=false, colormap="viridis")
            nexttile
            this.testObj.plot(r2, id="sms", combo=c2, use_cbrewer2=false, colormap="viridis")
            nexttile
            this.testObj.plot(r3, id="dan", combo=c3, use_cbrewer2=false, colormap="viridis")
            nexttile
            this.testObj.plot(r4, id="van", combo=c4, use_cbrewer2=false, colormap="viridis")
            nexttile
            this.testObj.plot(r5, id="limbic", combo=c5, use_cbrewer2=false, colormap="viridis")
            nexttile
            this.testObj.plot(r6, id="con", combo=c6, use_cbrewer2=false, colormap="viridis")
            nexttile
            this.testObj.plot(r7, id="dmn", combo=c7, use_cbrewer2=false, colormap="viridis")
            fontsize(scale=1.25);
        end
        function test_7rsns_mean_and_wm(this)
            [r1,c1] = this.testObj.parc("vis", statistic="mean");
            [r2,c2] = this.testObj.parc("sms", statistic="mean");
            [r3,c3] = this.testObj.parc("dan", statistic="mean");
            [r4,c4] = this.testObj.parc("van", statistic="mean");
            [r5,c5] = this.testObj.parc("limbic", statistic="mean");
            [r6,c6] = this.testObj.parc("con", statistic="mean");
            [r7,c7] = this.testObj.parc("dmn", statistic="mean");
            [rw,cw] = this.testObj.parc("white");

            cb = cbrewer2("Dark2", 7);

            tiledlayout(7, 1) %, TileSpacing="compact", Padding="compact")
            nexttile
            this.testObj.plot(rw, id="wm", use_cbrewer2=false, colormap="cividis", title="vis")
            hold on
            plot(r1, Color=cb(1,:), LineWidth=2)
            hold off
            nexttile
            this.testObj.plot(rw, id="wm", use_cbrewer2=false, colormap="cividis", title="sms")
            hold on
            plot(r2, Color=cb(2,:), LineWidth=2)
            hold off
            nexttile
            this.testObj.plot(rw, id="wm", use_cbrewer2=false, colormap="cividis", title="dan")
            hold on
            plot(r3, Color=cb(3,:), LineWidth=2)
            hold off
            nexttile
            this.testObj.plot(rw, id="wm", use_cbrewer2=false, colormap="cividis", title="van")
            hold on
            plot(r4, Color=cb(4,:), LineWidth=2)
            hold off
            nexttile
            this.testObj.plot(rw, id="wm", use_cbrewer2=false, colormap="cividis", title="limbic")
            hold on
            plot(r5, Color=cb(5,:), LineWidth=2)
            hold off
            nexttile
            this.testObj.plot(rw, id="wm", use_cbrewer2=false, colormap="cividis", title="con")
            hold on
            plot(r6, Color=cb(6,:), LineWidth=2)
            hold off
            nexttile
            this.testObj.plot(rw, id="wm", use_cbrewer2=false, colormap="cividis", title="dmn")
            hold on
            plot(r7, Color=cb(7,:), LineWidth=2)
            hold off
            fontsize(scale=1.25);
        end
    end
    
    methods (TestClassSetup)
        function setupTimeSeries(this)
            import mlvg.*
            cd("/Users/jjlee/Documents/PapersMine/PapersInProgress/PET_connectivity/CO");
            this.testObj_ = TimeSeries(this.pet, t0=60);
        end
    end
    
    methods (TestMethodSetup)
        function setupTimeSeriesTest(this)
            this.testObj = this.testObj_;
            this.addTeardown(@this.cleanTestMethod)
        end
    end
    
    properties (Access = private)
        testObj_
    end
    
    methods (Access = private)
        function cleanTestMethod(this)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
