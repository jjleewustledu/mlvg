runtests("mlvg_unittest.Test_Kudomi2018","ProcedureName","test_rho")
Running mlvg_unittest.Test_Kudomi2018
  Kudomi2018 with properties:

            Acond: []
            Bcond: []
             blur: 9
    diagonal_only: 0
               dt: 1
        Nsafeparc: 1276301
               Nx: 400
                p: 0.800000000000000
         safeparc: [176×208×176 mlfourd.ImagingContext2]
          scanner: [1×1 mlsiemens.BiographMMRDevice]
      sessionData: [1×1 mlraichle.SessionData]
           thresh: 0.100000000000000
            time0: 25
            timeF: 120
               Nt: 94
            times: [1×94 double]
         timesMid: [1×20 double]
              tol: 5.551115123125783e-15
      voxelVolume: 1.000000000000000e-03
          A_cache: [400×400 double]
          B_cache: []
    drho_dt_cache: [400×94 double]
        rho_cache: [400×94 double]
        tol_cache: []
          X_cache: []

75          end

            A = this.testObj.A();
            Ainv = pinv(A, this.testObj.tol);
            B = this.testObj.B();
            K1 = this.testObj.K1();
            
            figure; imagesc(A)
            figure; imagesc(Ainv)
            figure; imagesc(B)
            figure; plot(K1)