classdef Test_Kudomi2018 < matlab.unittest.TestCase
	%% TEST_KUDOMI2018 

	%  Usage:  >> results = run(mlvg_unittest.Test_Kudomi2018)
 	%          >> result  = run(mlvg_unittest.Test_Kudomi2018, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 13-Jul-2021 23:23:15 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlvg/test/+mlvg_unittest.
 	%% It was developed on Matlab 9.10.0.1669831 (R2021a) Update 2 for MACI64.  Copyright 2021 John Joowon Lee.
 	
	properties
        blur = 4.3
        diagonal_only = false
        dt = 1
        Nx = []
 		registry
        sessd
 		testObj
        time0 = 0
        timeF = 500
        tol = []
 	end

	methods (Test)
		function test_afun(this)
 			import mlvg.*;
 			this.assumeEqual(1,1);
 			this.verifyEqual(1,1);
 			this.assertEqual(1,1);
        end
        function test_aif(this)
            aifs = this.testObj.aifs();
            figure; surf(aifs + dipmax(aifs)/2);
            hold on
            imagesc(aifs)
            hold off
            aif = this.testObj.aif();
            figure; plot(mean(aif, 1))
            disp(this.testObj)
        end
        function test_ctor(this)
            disp(this.testObj)
        end
        function test_A_B(this)
            A = this.testObj.A();
            Ainv = pinv(A, this.testObj.tol);
            B = this.testObj.B();
            K1 = this.testObj.K1();
            
            figure; imagesc(A)
            figure; imagesc(Ainv)
            figure; imagesc(B)
            figure; plot(K1)
            disp(this.testObj)
            
            figure; plot(this.testObj.times, this.testObj.aifs())
        end
        function test_cond(this)
            this.testObj.A();            
            this.testObj.B();
            disp(this.testObj)
        end
        function test_K1(this)
            figure; plot(this.testObj.K1());
            disp(this.testObj)
        end
        function test_rho(this)
            rho = this.testObj.rho();
            drho_dt = this.testObj.drho_dt();
            d2rho_dt2 = this.testObj.d2rho_dt2();
            [Kvar,Kvar_] = this.testObj.Kvariation();
            
            figure; imagesc(rho)
            figure; imagesc(drho_dt)
            figure; imagesc(d2rho_dt2)
            figure; imagesc(Kvar)
            figure; imagesc(Kvar_)
            figure; plot(rho(50,:))
            figure; plot(drho_dt(50,:))
            figure; plot(d2rho_dt2(50,:))
            figure; plot(mean(Kvar,1))
            figure; plot(mean(Kvar_,1))
            disp(this.testObj)
        end
        function test_X(this)
            X = this.testObj.X();
            size(X)
            figure; imagesc(X)
            % scaxis([0 2])
        end
        function test_safeparc(this)
            this.testObj.safeparc.fsleyes
        end
	end

 	methods (TestClassSetup)
		function setupKudomi2018(this)
 			import mlvg.*;
            this.testObj_ = Kudomi2018( ...
                'sessionData', this.buildSessionData(), ...
                'blur', this.blur, ...
                'diagonal_only', this.diagonal_only, ...
                'dt', this.dt, ...
                'Nx', this.Nx, ...
                'time0', this.time0, ...
                'timeF', this.timeF, ...
                'tol', this.tol);
 		end
	end

 	methods (TestMethodSetup)
		function setupKudomi2018Test(this)
 			this.testObj = this.testObj_;
 			this.addTeardown(@this.cleanTestMethod);
 		end
	end

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
		function cleanTestMethod(this)
        end
        function sd = buildSessionData(this)
            import mlraichle.*;
            subd = SubjectData('subjectFolder', 'sub-S58163');
            sesfs = subd.subFolder2sesFolders('sub-S58163');
            sd = SessionData( ...
                'studyData', StudyData(), ...
                'projectData', ProjectData('sessionStr', sesfs{end}), ...
                'subjectData', subd, ...
                'sessionFolder', sesfs{end}, ...
                'scanIndex', 1, ...
                'tracer', 'HO', ...
                'ac', true, ...
                'region', 'wmparc1', ...
                'metric', 'fs');
        end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

