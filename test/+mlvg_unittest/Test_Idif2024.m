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

        function test_quotient(this)
            cd("/Volumes/PrecunealSSD/Singularity/CCIR_01211/derivatives")
            %fn_idif = "sub-all_ses-all_trc-ho_proc-schaefer-idif-raichleks-median.dtseries.nii";
            %fn_aif = "sub-all_ses-all_trc-ho_proc-schaefer-twilite-raichleks-median.dtseries.nii";
            %fn_idif = "sub-all_ses-all_trc-oo_proc-schaefer-idif-mintunks-median.dtseries.nii";
            %fn_aif = "sub-all_ses-all_trc-oo_proc-schaefer-twilite-mintunks-median.dtseries.nii";
            fn_idif = "sub-all_ses-all_trc-oo_proc-schaefer-idif-cmro2-median.dscalar.nii";
            fn_aif = "sub-all_ses-all_trc-oo_proc-schaefer-twilite-cmro2-median.dscalar.nii";
            c_idif = cifti_read(char(fn_idif));
            c_aif = cifti_read(char(fn_aif));
            c = c_idif;
            c.cdata = c_idif.cdata ./ c_aif.cdata;
            c.cdata(isnan(c.cdata)) = 0;
            cifti_write(c, char(strrep(fn_idif, "idif", "idif-over-aif")))
        end

        function test_martinv1(this)
            idif2024 = mlvg.Idif2024();
            pwd0 = pushd(idif2024.derivatives_path);

            % nii = idif2024.martinv1(sub="sub-108293");
            % imagesc(nii);
            % nii_idif = idif2024.martinv1();
            % imagesc(nii_idif);
            % nii_median = idif2024.martinv1(stats="mean");
            % imagesc(nii_median);
            % nii_iqr = idif2024.martinv1(stats="iqr");
            % imagesc(nii_iqr);
            % nii_coeffvar = nii_iqr ./ nii_median;
            % imagesc(nii_coeffvar);

            cii_idif = idif2024.martinv1(stats="mean", input_func="idif", typeclass="cifti");
            disp(cii_idif)
            cifti_write(cii_idif, convertStringsToChars(cii_idif.fqfn));

            cii_twil = idif2024.martinv1(stats="mean", input_func="twilite", typeclass="cifti");
            disp(cii_twil)
            cifti_write(cii_twil, convertStringsToChars(cii_twil.fqfn));

            cii_baresid = cii_idif;
            unbiased_est = 0.5*(cii_idif.cdata + cii_twil.cdata);
            cii_baresid.cdata = (cii_idif.cdata - cii_twil.cdata) ./ unbiased_est;
            cii_baresid.fqfn = strrep(cii_baresid.fqfn, "idif", "baresid");
            disp(cii_baresid)
            cifti_write(cii_baresid, convertStringsToChars(cii_baresid.fqfn))

            % write & view all subjects
            cii = idif2024.martinv1(typeclass="cifti");
            cifti_write(cii, convertStringsToChars(stackstr() + ".dtseries.nii"))
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
            % nii_median = idif2024.raichleks(stats="mean");
            % imagesc(nii_median);
            % nii_iqr = idif2024.raichleks(stats="iqr");
            % imagesc(nii_iqr);
            % nii_coeffvar = nii_iqr ./ nii_median;
            % imagesc(nii_coeffvar);

            cii_idif = idif2024.raichleks(stats="mean", input_func="idif", typeclass="cifti");
            disp(cii_idif)
            cifti_write(cii_idif, convertStringsToChars(cii_idif.fqfn));

            cii_twil = idif2024.raichleks(stats="mean", input_func="twilite", typeclass="cifti");
            disp(cii_twil)
            cifti_write(cii_twil, convertStringsToChars(cii_twil.fqfn));

            cii_baresid = cii_idif;
            unbiased_est = 0.5*(cii_idif.cdata + cii_twil.cdata);
            cii_baresid.cdata = (cii_idif.cdata - cii_twil.cdata) ./ unbiased_est;
            cii_baresid.fqfn = strrep(cii_baresid.fqfn, "idif", "baresid");
            disp(cii_baresid)
            cifti_write(cii_baresid, convertStringsToChars(cii_baresid.fqfn))

            % write & view all subjects
            cii = idif2024.raichleks(typeclass="cifti");
            cifti_write(cii, convertStringsToChars(stackstr() + ".dtseries.nii"))
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
            % nii_median = idif2024.mintunks(stats="mean");
            % imagesc(nii_median);
            % nii_iqr = idif2024.mintunks(stats="iqr");
            % imagesc(nii_iqr);
            % nii_coeffvar = nii_iqr ./ nii_median;
            % imagesc(nii_coeffvar);

            cii_idif = idif2024.mintunks(stats="mean", input_func="idif", typeclass="cifti");
            disp(cii_idif)
            cifti_write(cii_idif, convertStringsToChars(cii_idif.fqfn));

            cii_twil = idif2024.mintunks(stats="mean", input_func="twilite", typeclass="cifti");
            disp(cii_twil)
            cifti_write(cii_twil, convertStringsToChars(cii_twil.fqfn));

            cii_baresid = cii_idif;
            unbiased_est = 0.5*(cii_idif.cdata + cii_twil.cdata);
            cii_baresid.cdata = (cii_idif.cdata - cii_twil.cdata) ./ unbiased_est;
            cii_baresid.fqfn = strrep(cii_baresid.fqfn, "idif", "baresid");
            disp(cii_baresid)
            cifti_write(cii_baresid, convertStringsToChars(cii_baresid.fqfn))

            % write & view all subjects
            cii = idif2024.mintunks(typeclass="cifti");
            cifti_write(cii, convertStringsToChars(stackstr() + ".dtseries.nii"))
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
            % nii_median = idif2024.huangks(stats="mean");
            % imagesc(nii_median);
            % nii_iqr = idif2024.huangks(stats="iqr");
            % imagesc(nii_iqr);
            % nii_coeffvar = nii_iqr ./ nii_median;
            % imagesc(nii_coeffvar);

            cii_idif = idif2024.huangks(stats="mean", input_func="idif", typeclass="cifti");
            disp(cii_idif)
            cifti_write(cii_idif, convertStringsToChars(cii_idif.fqfn));

            cii_twil = idif2024.huangks(stats="mean", input_func="twilite", typeclass="cifti");
            disp(cii_twil)
            cifti_write(cii_twil, convertStringsToChars(cii_twil.fqfn));

            cii_baresid = cii_idif;
            unbiased_est = 0.5*(cii_idif.cdata + cii_twil.cdata);
            cii_baresid.cdata = (cii_idif.cdata - cii_twil.cdata) ./ unbiased_est;
            cii_baresid.fqfn = strrep(cii_baresid.fqfn, "idif", "baresid");
            disp(cii_baresid)
            cifti_write(cii_baresid, convertStringsToChars(cii_baresid.fqfn))

            % write & view all subjects
            cii = idif2024.huangks(typeclass="cifti");
            cifti_write(cii, convertStringsToChars(stackstr() + ".dtseries.nii"))
            % mysystem("wb_view")

            popd(pwd0);
        end

        function test_cmro2(this)
            idif2024 = mlvg.Idif2024();
            pwd0 = pushd(idif2024.derivatives_path);

            % nii = idif2024.cmro2(stats="", input_func="idif", typeclass="nifti");
            % imagesc(nii);
            
            cii_idif = idif2024.cmro2(stats="mean", input_func="idif", typeclass="cifti");
            disp(cii_idif)
            cifti_write(cii_idif, convertStringsToChars(cii_idif.fqfn));

            cii_twil = idif2024.cmro2(stats="mean", input_func="twilite", typeclass="cifti");
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

        function test_cmrglc(this)
            idif2024 = mlvg.Idif2024();
            pwd0 = pushd(idif2024.derivatives_path);

            nii = idif2024.cmrglc(stats="", input_func="idif", typeclass="nifti");
            imagesc(nii);
            
            cii_idif = idif2024.cmrglc(stats="mean", input_func="idif", typeclass="cifti");
            disp(cii_idif)
            cifti_write(cii_idif, convertStringsToChars(cii_idif.fqfn));

            cii_twil = idif2024.cmrglc(stats="mean", input_func="twilite", typeclass="cifti");
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

            % return

            idif2024 = mlvg.Idif2024();
            pwd0 = pushd(idif2024.derivatives_path);

            %model = "Raichle1983";
            model = "Mintun1984";
            %model = "Huang1980";
            %trc = "trc-ho";
            trc = "trc-oo";
            %trc = "trc-fdg";
            sub = "sub-*";
            centrality = "median";
            deviation = "iqr";

            nii_idif = idif2024.logz(stats=centrality, sub=sub, input_func="idif", model=model, trc=trc);
            disp(nii_idif);

            nii_twil = idif2024.logz(stats=centrality, sub=sub, input_func="twil", model=model, trc=trc);
            disp(nii_twil);

            %% median, iqr, coeff. var.

            nii_centrality = idif2024.logz(stats=centrality, input_func="idif", model=model, trc=trc);
            plot(nii_centrality);

            nii_deviation = idif2024.logz(stats=deviation, input_func="idif", model=model, trc=trc);
            plot(nii_deviation);

            nii_coeffvar = nii_deviation ./ nii_centrality;
            plot(nii_coeffvar);

            %% log(Bayes factor) ~ Dlogz

            nii_bf = nii_idif - nii_twil;
            plot(nii_bf);

            %% cifti to write

            cii_idif = idif2024.logz(stats=centrality, input_func="idif", typeclass="cifti", model=model, trc=trc);
            disp(cii_idif)
            cifti_write(cii_idif, convertStringsToChars(cii_idif.fqfn));

            cii_twil = idif2024.logz(stats=centrality, input_func="twilite", typeclass="cifti", model=model, trc=trc);
            disp(cii_twil)
            cifti_write(cii_twil, convertStringsToChars(cii_twil.fqfn));

            cii_bf = idif2024.as(nii_bf, "cifti");
            s = split(cii_idif.fqfn, "_proc-");
            cii_bf.fqfn = s(1) + "_proc-schaefer-bfactors-" + centrality + ".dscalar.nii";
            disp(cii_bf)
            cifti_write(cii_bf, convertStringsToChars(cii_bf.fqfn));

            popd(pwd0);
        end

        function test_bayes_factors(this)

            return

            idif2024 = mlvg.Idif2024();
            pwd0 = pushd(idif2024.derivatives_path);

            %model = "Raichle1983"; trc = "trc-ho";
            model = "Mintun1984"; trc = "trc-oo";
            %model = "Huang1980"; trc = "trc-fdg";           

            nii_bf = idif2024.bayes_factors( ...
                typeclass="nifti", stats="mean", sub="sub-108293", model=model, trc=trc);
            disp(nii_bf);
            imagesc(nii_bf);

            cii_bf = idif2024.as(nii_bf, "cifti");
            s = split(cii_bf.fqfn, "_proc-");
            cii_bf.fqfn = s(1) + "_proc-schaefer-" + model + "-bfactors-mean.dscalar.nii";
            disp(cii_bf)
            cifti_write(cii_bf, convertStringsToChars(cii_bf.fqfn));

            popd(pwd0);
        end
        function test_rm_raincloud(this)

            return

            idif2024 = mlvg.Idif2024();

            pwd0 = pushd(idif2024.derivatives_path);

            metrics = [
                "OEF", ...
                "CMRO2", ...
                "OGI"];
            axis_labels = [...
                "OEF", ...
                "CMRO2 (mmol L^{-1} min^{-1})", ...
                "OGI"];
            for idx = 1:numel(metrics)
                idif2024.build_rm_raincloud( ...
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
