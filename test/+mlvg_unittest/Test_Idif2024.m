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

            nii = idif2024.martinv1(sub="sub-108293");
            imagesc(nii);
            nii_idif = idif2024.martinv1();
            imagesc(nii_idif);
            nii_median = idif2024.martinv1(stats="median");
            imagesc(nii_median);
            nii_iqr = idif2024.martinv1(stats="iqr");
            imagesc(nii_iqr);
            nii_coeffvar = nii_iqr ./ nii_median;
            imagesc(nii_coeffvar);

            cii_idif = idif2024.martinv1(stats="median", input_func="idif", typeclass="cifti");
            disp(cii_idif)
            %cifti_write(cii_idif, convertStringsToChars(cii_idif.fqfn));

            cii_twil = idif2024.martinv1(stats="median", input_func="twilite", typeclass="cifti");
            disp(cii_twil)
            %cifti_write(cii_twil, convertStringsToChars(cii_twil.fqfn));

            cii_didif = cii_idif;
            cii_didif.cdata = cii_idif.cdata - cii_twil.cdata;
            cii_didif.fqfn = strrep(cii_didif.fqfn, "idif", "didif");
            disp(cii_didif)
            %cifti_write(cii_didif, convertStringsToChars(cii_didif.fqfn))

            % write & view all subjects
            %cii = idif2024.martinv1(typeclass="cifti");
            %cifti_write(cii, convertStringsToChars(stackstr() + ".dtseries.nii"))
            %mysystem("wb_view")

            popd(pwd0);
        end
        function test_raichleks(this)
            idif2024 = mlvg.Idif2024();
            pwd0 = pushd(idif2024.derivatives_path);

            nii = idif2024.raichleks(sub="sub-108293");
            imagesc(nii);
            nii_idif = idif2024.raichleks();
            view(nii_idif);
            nii_median = idif2024.raichleks(stats="median");
            imagesc(nii_median);
            nii_iqr = idif2024.raichleks(stats="iqr");
            imagesc(nii_iqr);
            nii_coeffvar = nii_iqr ./ nii_median;
            imagesc(nii_coeffvar);

            cii_idif = idif2024.raichleks(stats="median", input_func="idif", typeclass="cifti");
            disp(cii_idif)
            %cifti_write(cii_idif, convertStringsToChars(cii_idif.fqfn));

            cii_twil = idif2024.raichleks(stats="median", input_func="twilite", typeclass="cifti");
            disp(cii_twil)
            %cifti_write(cii_twil, convertStringsToChars(cii_twil.fqfn));

            cii_didif = cii_idif;
            cii_didif.cdata = cii_idif.cdata - cii_twil.cdata;
            cii_didif.fqfn = strrep(cii_didif.fqfn, "idif", "didif");
            disp(cii_didif)
            %cifti_write(cii_didif, convertStringsToChars(cii_didif.fqfn))

            % write & view all subjects
            %cii = idif2024.raichleks(typeclass="cifti");
            %cifti_write(cii, convertStringsToChars(stackstr() + ".dtseries.nii"))
            %mysystem("wb_view")

            popd(pwd0);
        end
        function test_mintunks(this)
            idif2024 = mlvg.Idif2024();
            pwd0 = pushd(idif2024.derivatives_path);

            nii = idif2024.mintunks(sub="sub-108293");
            imagesc(nii);
            nii_idif = idif2024.mintunks();
            view(nii_idif);
            nii_median = idif2024.mintunks(stats="median");
            imagesc(nii_median);
            nii_iqr = idif2024.mintunks(stats="iqr");
            imagesc(nii_iqr);
            nii_coeffvar = nii_iqr ./ nii_median;
            imagesc(nii_coeffvar);

            cii_idif = idif2024.mintunks(stats="median", input_func="idif", typeclass="cifti");
            disp(cii_idif)
            cifti_write(cii_idif, convertStringsToChars(cii_idif.fqfn));

            cii_twil = idif2024.mintunks(stats="median", input_func="twilite", typeclass="cifti");
            disp(cii_twil)
            cifti_write(cii_twil, convertStringsToChars(cii_twil.fqfn));

            cii_didif = cii_idif;
            cii_didif.cdata = cii_idif.cdata - cii_twil.cdata;
            cii_didif.fqfn = strrep(cii_didif.fqfn, "idif", "didif");
            disp(cii_didif)
            cifti_write(cii_didif, convertStringsToChars(cii_didif.fqfn))

            % write & view all subjects
            %cii = idif2024.mintunks(typeclass="cifti");
            %cifti_write(cii, convertStringsToChars(stackstr() + ".dtseries.nii"))
            %mysystem("wb_view")
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
