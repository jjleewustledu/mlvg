classdef Test_PowerSpectralDensities < matlab.unittest.TestCase
    %% https://claude.ai/chat/c44fe792-1083-48d9-a879-e2757bfa0236
    %  
    %  Created 10-Jul-2024 12:35:47 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/test/+mlvg_unittest.
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

        function test_cwt(this)
            % Load your data
            data = this.testObj.grey;
            t = this.testObj.adj_times;
            samplingf = this.testObj.sampling_freq;

            % Perform CWT
            [wt, f] = cwt(data, 'amor', samplingf);

            % Plot scalogram
            figure;
            pcolor(t, f, abs(wt).^2);
            shading flat;
            set(gca, 'YScale', 'log');
            xlabel('Time');
            ylabel('Frequency');
            title('Wavelet Scalogram');
            colorbar;
        end

        function test_gws(this)
            % Load your data
            data = this.testObj.grey;
            t = this.testObj.adj_times;
            samplingf = this.testObj.sampling_freq;

            % Perform CWT
            [wt, f] = cwt(data, 'amor', samplingf);

            % Compute global wavelet spectrum
            gws = mean(abs(wt).^2, 2);

            % Plot global wavelet spectrum
            figure;
            loglog(f, gws);
            xlabel('Frequency');
            ylabel('Power');
            title('Global Wavelet Spectrum');
        end

        function test_powerlaws(this)
            % Load your data
            data = this.testObj.grey;
            t = this.testObj.adj_times;
            samplingf = this.testObj.sampling_freq;

            % Perform CWT
            [wt, f] = cwt(data, 'amor', samplingf);

            % Compute global wavelet spectrum
            gws = mean(abs(wt).^2, 2);

            % Identify linear regions in log-log plot
            [p, S] = polyfit(log10(f), log10(gws), 1);
            slope = p(1);

            hold on;
            loglog(f, 10.^polyval(p, log10(f)), 'r--');
            legend('Global Wavelet Spectrum', 'Fitted Line');
            text(0.1, max(gws), sprintf('Slope: %.2f', slope));
        end

        function test_wavelet_coherence(this)
        end

        function test_dwt(this)
            % Perform DWT
            data = this.testObj.grey;
            wname = 'db4'; % Choose an appropriate wavelet
            level = 5; % Choose an appropriate decomposition level
            [c, l] = wavedec(data, level, wname);

            % Compute and plot approximation and detail coefficients
            cd5 = c(l(1)+1:l(1)+l(2)); %#ok<*NASGU>
            cd4 = c(l(1)+l(2)+1:l(1)+l(2)+l(3));
            cd3 = c(l(1)+l(2)+l(3)+1:l(1)+l(2)+l(3)+l(4));
            cd2 = c(l(1)+l(2)+l(3)+l(4)+1:l(1)+l(2)+l(3)+l(4)+l(5));
            cd1 = c(l(1)+l(2)+l(3)+l(4)+l(5)+1:l(1)+l(2)+l(3)+l(4)+l(5)+l(6));

            figure;
            subplot(level+1, 1, 1);
            plot(data);
            title('Original Signal');

            for i = 1:level
                subplot(level+1, 1, i+1);
                eval(['plot(cd', num2str(level+1-i), ');']);
                title(['Detail Coefficients - Level ', num2str(level+1-i)]);
            end

            % Compute energy at each level
            energies = zeros(1, level+1);
            for i = 1:level
                d = detcoef(c, l, i);
                energies(i) = sum(d.^2);
            end
            energies(level+1) = sum(appcoef(c, l, wname).^2);

            % Plot energy distribution
            figure;
            bar(energies);
            xlabel('Decomposition Level');
            ylabel('Energy');
            title('Energy Distribution Across Wavelet Scales');
        end
    end
    
    methods (TestClassSetup)
        function setupPowerSpectralDensities(this)
            import mlvg.*
            cd("/Users/jjlee/Documents/PapersMine/PapersInProgress/PET_connectivity/CO");
            this.testObj_ = PowerSpectralDensities(this.pet, t0=60);
        end
    end
    
    methods (TestMethodSetup)
        function setupPowerSpectralDensitiesTest(this)
            this.testObj = copy(this.testObj_);
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
