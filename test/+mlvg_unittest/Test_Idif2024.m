classdef Test_Idif2024 < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 06-Mar-2024 23:29:15 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/test/+mlvg_unittest.
    %  Developed on Matlab 23.2.0.2515942 (R2023b) Update 7 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mlvg.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end

        function test_martinv1(this)
            idif2024 = mlvg.Idif2024();
            pwd0 = pushd(idif2024.derivatives_path);

            % nii = idif2024.martinv1(sub="sub-108293");
            % imagesc(nii);
            % nii_idif = idif2024.martinv1();
            % imagesc(nii_idif);
            % nii_median = idif2024.martinv1(stats="median");
            % imagesc(nii_median);
            % nii_iqr = idif2024.martinv1(stats="iqr");
            % imagesc(nii_iqr);
            % nii_coeffvar = nii_iqr ./ nii_median;
            % imagesc(nii_coeffvar);

            cii_idif = idif2024.martinv1(stats="median", input_func="idif", typeclass="cifti");
            disp(cii_idif)
            % cifti_write(cii_idif, convertStringsToChars(cii_idif.fqfn));

            cii_twil = idif2024.martinv1(stats="median", input_func="twilite", typeclass="cifti");
            disp(cii_twil)
            % cifti_write(cii_twil, convertStringsToChars(cii_twil.fqfn));

            cii_baresid = cii_idif;
            unbiased_est = 0.5*(cii_idif.cdata + cii_twil.cdata);
            cii_baresid.cdata = (cii_idif.cdata - cii_twil.cdata) ./ unbiased_est;
            cii_baresid.fqfn = strrep(cii_baresid.fqfn, "idif", "baresid");
            disp(cii_baresid)
            cifti_write(cii_baresid, convertStringsToChars(cii_baresid.fqfn))

            % write & view all subjects
            % cii = idif2024.martinv1(typeclass="cifti");
            % cifti_write(cii, convertStringsToChars(stackstr() + ".dtseries.nii"))
            % mysystem("wb_view")

            popd(pwd0);
        end

        function test_raichleks(this)
            idif2024 = mlvg.Idif2024();
            pwd0 = pushd(idif2024.derivatives_path);

            % nii = idif2024.raichleks(sub="sub-108293");
            % imagesc(nii);
            % nii_idif = idif2024.raichleks();
            % view(nii_idif);
            % nii_median = idif2024.raichleks(stats="median");
            % imagesc(nii_median);
            % nii_iqr = idif2024.raichleks(stats="iqr");
            % imagesc(nii_iqr);
            % nii_coeffvar = nii_iqr ./ nii_median;
            % imagesc(nii_coeffvar);

            cii_idif = idif2024.raichleks(stats="median", input_func="idif", typeclass="cifti");
            disp(cii_idif)
            % cifti_write(cii_idif, convertStringsToChars(cii_idif.fqfn));

            cii_twil = idif2024.raichleks(stats="median", input_func="twilite", typeclass="cifti");
            disp(cii_twil)
            % cifti_write(cii_twil, convertStringsToChars(cii_twil.fqfn));

            cii_baresid = cii_idif;
            unbiased_est = 0.5*(cii_idif.cdata + cii_twil.cdata);
            cii_baresid.cdata = (cii_idif.cdata - cii_twil.cdata) ./ unbiased_est;
            cii_baresid.fqfn = strrep(cii_baresid.fqfn, "idif", "baresid");
            disp(cii_baresid)
            cifti_write(cii_baresid, convertStringsToChars(cii_baresid.fqfn))

            % write & view all subjects
            % cii = idif2024.raichleks(typeclass="cifti");
            % cifti_write(cii, convertStringsToChars(stackstr() + ".dtseries.nii"))
            % mysystem("wb_view")

            popd(pwd0);
        end

        function test_mintunks(this)
            idif2024 = mlvg.Idif2024();
            pwd0 = pushd(idif2024.derivatives_path);

            % nii = idif2024.mintunks(sub="sub-108293");
            % imagesc(nii);
            % nii_idif = idif2024.mintunks();
            % view(nii_idif);
            % nii_median = idif2024.mintunks(stats="median");
            % imagesc(nii_median);
            % nii_iqr = idif2024.mintunks(stats="iqr");
            % imagesc(nii_iqr);
            % nii_coeffvar = nii_iqr ./ nii_median;
            % imagesc(nii_coeffvar);

            cii_idif = idif2024.mintunks(stats="median", input_func="idif", typeclass="cifti");
            disp(cii_idif)
            % cifti_write(cii_idif, convertStringsToChars(cii_idif.fqfn));

            cii_twil = idif2024.mintunks(stats="median", input_func="twilite", typeclass="cifti");
            disp(cii_twil)
             % cifti_write(cii_twil, convertStringsToChars(cii_twil.fqfn));

            cii_baresid = cii_idif;
            unbiased_est = 0.5*(cii_idif.cdata + cii_twil.cdata);
            cii_baresid.cdata = (cii_idif.cdata - cii_twil.cdata) ./ unbiased_est;
            cii_baresid.fqfn = strrep(cii_baresid.fqfn, "idif", "baresid");
            disp(cii_baresid)
            cifti_write(cii_baresid, convertStringsToChars(cii_baresid.fqfn))

            % write & view all subjects
            % cii = idif2024.mintunks(typeclass="cifti");
            % cifti_write(cii, convertStringsToChars(stackstr() + ".dtseries.nii"))
            % mysystem("wb_view")

            popd(pwd0);
        end

        function test_huangks(this)
            idif2024 = mlvg.Idif2024();
            pwd0 = pushd(idif2024.derivatives_path);

            % nii = idif2024.huangks(sub="sub-108293");
            % imagesc(nii);
            % nii_idif = idif2024.huangks();
            % view(nii_idif);
            % nii_median = idif2024.huangks(stats="median");
            % imagesc(nii_median);
            % nii_iqr = idif2024.huangks(stats="iqr");
            % imagesc(nii_iqr);
            % nii_coeffvar = nii_iqr ./ nii_median;
            % imagesc(nii_coeffvar);

            cii_idif = idif2024.huangks(stats="median", input_func="idif", typeclass="cifti");
            disp(cii_idif)
            % cifti_write(cii_idif, convertStringsToChars(cii_idif.fqfn));

            cii_twil = idif2024.huangks(stats="median", input_func="twilite", typeclass="cifti");
            disp(cii_twil)
            % cifti_write(cii_twil, convertStringsToChars(cii_twil.fqfn));

            cii_baresid = cii_idif;
            unbiased_est = 0.5*(cii_idif.cdata + cii_twil.cdata);
            cii_baresid.cdata = (cii_idif.cdata - cii_twil.cdata) ./ unbiased_est;
            cii_baresid.fqfn = strrep(cii_baresid.fqfn, "idif", "baresid");
            disp(cii_baresid)
            cifti_write(cii_baresid, convertStringsToChars(cii_baresid.fqfn))

            % write & view all subjects
            % cii = idif2024.huangks(typeclass="cifti");
            % cifti_write(cii, convertStringsToChars(stackstr() + ".dtseries.nii"))
            % mysystem("wb_view")

            popd(pwd0);
        end

        function test_cmro2(this)
            idif2024 = mlvg.Idif2024();
            pwd0 = pushd(idif2024.derivatives_path);

            % nii = idif2024.cmro2(stats="", input_func="idif", typeclass="nifti");
            % imagesc(nii);
            
            cii_idif = idif2024.cmro2(stats="median", input_func="idif", typeclass="cifti");
            disp(cii_idif)
            % cifti_write(cii_idif, convertStringsToChars(cii_idif.fqfn));

            cii_twil = idif2024.cmro2(stats="median", input_func="twilite", typeclass="cifti");
            disp(cii_twil)
            % cifti_write(cii_twil, convertStringsToChars(cii_twil.fqfn));

            cii_baresid = cii_idif;
            unbiased_est = 0.5*(cii_idif.cdata + cii_twil.cdata);
            cii_baresid.cdata = (cii_idif.cdata - cii_twil.cdata) ./ unbiased_est;
            cii_baresid.fqfn = strrep(cii_baresid.fqfn, "idif", "baresid");
            disp(cii_baresid)
            cifti_write(cii_baresid, convertStringsToChars(cii_baresid.fqfn))

            popd(pwd0);
        end

        function test_cmrglc(this)
            idif2024 = mlvg.Idif2024();
            pwd0 = pushd(idif2024.derivatives_path);

            nii = idif2024.cmrglc(stats="", input_func="idif", typeclass="nifti");
            imagesc(nii);
            
            cii_idif = idif2024.cmrglc(stats="median", input_func="idif", typeclass="cifti");
            disp(cii_idif)
            cifti_write(cii_idif, convertStringsToChars(cii_idif.fqfn));

            cii_twil = idif2024.cmrglc(stats="median", input_func="twilite", typeclass="cifti");
            disp(cii_twil)
            cifti_write(cii_twil, convertStringsToChars(cii_twil.fqfn));

            cii_baresid = cii_idif;
            unbiased_est = 0.5*(cii_idif.cdata + cii_twil.cdata);
            cii_baresid.cdata = (cii_idif.cdata - cii_twil.cdata) ./ unbiased_est;
            cii_baresid.fqfn = strrep(cii_baresid.fqfn, "idif", "baresid");
            disp(cii_baresid)
            cifti_write(cii_baresid, convertStringsToChars(cii_baresid.fqfn))

            popd(pwd0);
        end

        function test_ogi(this)
            idif2024 = mlvg.Idif2024();
            pwd0 = pushd(idif2024.derivatives_path);

            nii = idif2024.ogi(input_func="idif", typeclass="nifti");
            imagesc(nii);
            
            cii_idif = idif2024.ogi(input_func="idif", typeclass="cifti");
            disp(cii_idif)
            cifti_write(cii_idif, convertStringsToChars(cii_idif.fqfn));

            cii_twil = idif2024.ogi(input_func="twilite", typeclass="cifti");
            disp(cii_twil)
            cifti_write(cii_twil, convertStringsToChars(cii_twil.fqfn));

            cii_baresid = cii_idif;
            unbiased_est = 0.5*(cii_idif.cdata + cii_twil.cdata);
            cii_baresid.cdata = (cii_idif.cdata - cii_twil.cdata) ./ unbiased_est;
            cii_baresid.fqfn = strrep(cii_baresid.fqfn, "idif", "baresid");
            disp(cii_baresid)
            cifti_write(cii_baresid, convertStringsToChars(cii_baresid.fqfn))

            popd(pwd0);
        end

        function test_agi(this)
            idif2024 = mlvg.Idif2024();
            pwd0 = pushd(idif2024.derivatives_path);

            nii = idif2024.agi(input_func="idif", typeclass="nifti");
            imagesc(nii);
            
            cii_idif = idif2024.agi(input_func="idif", typeclass="cifti");
            disp(cii_idif)
            cifti_write(cii_idif, convertStringsToChars(cii_idif.fqfn));

            cii_twil = idif2024.agi(input_func="twilite", typeclass="cifti");
            disp(cii_twil)
            cifti_write(cii_twil, convertStringsToChars(cii_twil.fqfn));

            cii_baresid = cii_idif;
            unbiased_est = 0.5*(cii_idif.cdata + cii_twil.cdata);
            cii_baresid.cdata = (cii_idif.cdata - cii_twil.cdata) ./ unbiased_est;
            cii_baresid.fqfn = strrep(cii_baresid.fqfn, "idif", "baresid");
            disp(cii_baresid)
            cifti_write(cii_baresid, convertStringsToChars(cii_baresid.fqfn))

            popd(pwd0);
        end

        function test_logz(this)
            idif2024 = mlvg.Idif2024();
            pwd0 = pushd(idif2024.derivatives_path);

            %model = "Raichle1983";
            model = "Mintun1984";
            %model = "Huang1980";
            %trc = "trc-ho";
            trc = "trc-oo";
            %trc = "trc-fdg";

            nii_idif = idif2024.logz(stats="median", sub="sub-108293", input_func="idif", model=model, trc=trc);
            disp(nii_idif);

            nii_twil = idif2024.logz(stats="median", sub="sub-108293", input_func="twil", model=model, trc=trc);
            disp(nii_twil);

            %% median, iqr, coeff. var.

            nii_median = idif2024.logz(stats="median", input_func="idif", model=model, trc=trc);
            plot(nii_median);

            nii_iqr = idif2024.logz(stats="iqr", input_func="idif", model=model, trc=trc);
            plot(nii_iqr);

            nii_coeffvar = nii_iqr ./ nii_median;
            plot(nii_coeffvar);

            %% log(Bayes factor) ~ Dlogz

            nii_bf = nii_idif - nii_twil;
            plot(nii_bf);

            %% cifti to write

            cii_idif = idif2024.logz(stats="median", input_func="idif", typeclass="cifti", model=model, trc=trc);
            disp(cii_idif)
            cifti_write(cii_idif, convertStringsToChars(cii_idif.fqfn));

            cii_twil = idif2024.logz(stats="median", input_func="twilite", typeclass="cifti", model=model, trc=trc);
            disp(cii_twil)
            cifti_write(cii_twil, convertStringsToChars(cii_twil.fqfn));

            cii_bf = idif2024.as(nii_bf, "cifti");
            s = split(cii_idif.fqfn, "_proc-");
            cii_bf.fqfn = s(1) + "_proc-schaefer-bfactors-median.dscalar.nii";
            disp(cii_bf)
            cifti_write(cii_bf, convertStringsToChars(cii_bf.fqfn));

            popd(pwd0);
        end

        function test_bayes_factors(this)
            idif2024 = mlvg.Idif2024();
            pwd0 = pushd(idif2024.derivatives_path);

            model = "Raichle1983"; trc = "trc-ho";
            %model = "Mintun1984"; trc = "trc-oo";
            %model = "Huang1980"; trc = "trc-fdg";           

            nii_bf = idif2024.bayes_factors( ...
                typeclass="nifti", stats="median", sub="sub-108293", model=model, trc=trc);
            disp(nii_bf);
            imagesc(nii_bf);

            cii_bf = idif2024.as(nii_bf, "cifti");
            s = split(cii_bf.fqfn, "_proc-");
            cii_bf.fqfn = s(1) + "_proc-schaefer-" + model + "-bfactors-median.dscalar.nii";
            disp(cii_bf)
            cifti_write(cii_bf, convertStringsToChars(cii_bf.fqfn));

            popd(pwd0);
        end
        function test_rm_raincloud(this)
            idif2024 = mlvg.Idif2024();
            pwd0 = pushd(idif2024.derivatives_path);

            metrics = [ ...
                "CBV", ...
                "CBF", ...
                "lambda", ...
                "PS", ...
                "OEF", ...
                "vcapillary", ...
                "fwatermetab", ...
                "CMRO2", ...
                "upper_k1", ...
                "k2", ...
                "k3", ...
                "k4", ...
                "CMRglc", ...
                "OGI", ...
                "AGI", ...
                "logZ_ho", ...
                "logZ_oo", ...
                "logZ_fdg"];
            axis_labels = [...
                "CBV (mL cm^{-3})", ...
                "CBF (mL cm^{-3} min^{-1})", ...
                "\lambda (mL cm^{-3})", ...
                "PS (mL cm^{-3} min^{-1})", ...
                "OEF", ...
                "V_{post,cap}", ...
                "f_{water}", ...
                "CMRO2 (mmol L^{-1} min^{-1})", ... 
                "K_1 (mL cm^{-3} min^{-1})", ...
                "k_2 (min^{-1})", ...
                "k_3 (min^{-1})", ...
                "k_4 (min^{-1})", ...
                "CMRglc (mmol L^{-1} min^{-1})", ...
                "OGI", ...
                "AGI (mmol L^{-1} min^{-1})", ...
                "log(Z) [^{15}O] H_2O", ...
                "log(Z) [^{15}O] O_2", ...
                "log(Z) [^{18}F] FDG"];
            for idx = 9:numel(metrics)
                idif2024.rm_raincloud( ...
                    metric=metrics(idx), axis_label=axis_labels(idx));
            end

            popd(pwd0);
        end
    end
    
    methods (TestClassSetup)
        function setupIdif2024(this)
            import mlvg.*
            this.testObj_ = Idif2024();
        end
    end
    
    methods (TestMethodSetup)
        function setupIdif2024Test(this)
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
