classdef Clinical < handle
    %% line1
    %  line2
    %  
    %  Created 06-Sep-2025 15:37:21 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.
    
    methods
        function this = Clinical(csv_fqfn)
            arguments
                csv_fqfn {mustBeFile} = ...
                    fullfile(getenv("HOME"), "MATLAB-Drive", "mlvg", "data", ...
                    "VGLabAgingDementia-VisionScanDataForVla_DATA_2025-09-06_1513.csv")
            end
            
            T = readtable(csv_fqfn);
            
            % simplify variables
            T.record_id = uint64(T.record_id);
            T.redcap_event_name = string(T.redcap_event_name);
            T.tau_scan_id = string(T.tau_scan_id);

            selection = ~isnan(T.sex_v2);
            T.sex(selection) = T.sex_v2(selection);  % male == 1; female == 0
            selection = ~isnan(T.sex_mag);
            T.sex(selection) = T.sex_mag(selection);
            selection = ~isnan(T.race_v2);
            T.race(selection) = T.race_v2(selection);
            selection = ~isnan(T.race_mag);
            T.race(selection) = T.race_mag(selection);
            selection = ~isnan(T.ethnicity_v2);
            T.ethnicity(selection) = T.ethnicity_v2(selection);
            selection = ~isnan(T.ethnicity_mag);
            T.ethnicity(selection) = T.ethnicity_mag(selection);

            pseudo_dob = T.pet_ag_date - years(T.age_pet_ag);
            T = addvars(T, pseudo_dob);
            pseudo_dob1 = T.pet_ag_date_wmh - years(T.age_pet_ag_wmh);
            selection = ~isnat(pseudo_dob1);
            T.pseudo_dob(selection) = pseudo_dob1(selection);
            pseudo_dob1 = T.date_of_session - years(T.age_pet_mag);
            selection = ~isnat(pseudo_dob1);
            T.pseudo_dob(selection) = pseudo_dob1(selection);
            T = renamevars(T, "date_of_session", "pet_mag_date");

            pet_metab_last_date = max([T.pet_ag_date, T.pet_ag_date_wmh, T.pet_mag_date], [], 2);
            T = addvars(T, pet_metab_last_date);

            T = removevars(T, ["sex_v2", "sex_mag", "race_v2", "race_mag", "ethnicity_v2", "ethnicity_mag", ...
                "age_pet_ag", "age_pet_ag_wmh", "age_pet_mag"]);

            this.table_internal_ = T;
        end

        function T = impute_hgb(~, T)
            missing = isnan(T.hgb);
            ismale = 1 == T.sex;
            isfem = 0 == T.sex;
            T.hgb(missing & ismale) = mlvg.Laboratory.Nijboer_2007("m");
            T.hgb(missing & isfem) = mlvg.Laboratory.Nijboer_2007("f");
        end

        function g = sub_to_age(this, sub, pet_date)
            %% returns age at PET scan

            if ~isnumeric(sub)
                sub = double(extractAfter(sub, "sub-"));
            end
            assert(isdatetime(pet_date))

            T = table(this);
            dob = T.pseudo_dob(sub == T.record_id, :);
            dob = mean(dob);
            g = years(pet_date - dob);
        end

        function g = sub_to_glc(this, sub)
            %% returns \mu mol/mL, most recent

            if ~isnumeric(sub)
                sub = double(extractAfter(sub, "sub-"));
            end

            Tlab = table_with_labs(this);
            glc = Tlab.glc_mg_dL(Tlab.record_id == sub & ~isnan(Tlab.glc_mg_dL));  % mg/dL
            if ~isempty(glc)
                g = this.molar_density_glc(glc(end));  % \mu mol/mL
            else 
                g = nan;
            end
        end

        function g = sub_to_o2_content(this, sub)
            %% returns \mu mol/mL, most recent

            if ~isnumeric(sub)
                sub = double(extractAfter(sub, "sub-"));
            end

            Tlab = table_with_labs(this);
            SpO2 = Tlab.PulseOx(Tlab.record_id == sub & ~isnan(Tlab.PulseOx));  % SpO2 \in [0, 100]
            hgb = Tlab.hgb(Tlab.record_id == sub & ~isnan(Tlab.hgb));  % hgb \in [5, 25]
            if ~isempty(SpO2) && ~isempty(hgb)
                CaO2 = this.Marino_eq_10_4(hgb(end), SpO2(end) * 0.01);  % mL/dL
                g = this.molar_density_O2(CaO2);
            else
                g = nan;
            end
        end

        function T = table(this, varargin)
            T = this.table_internal_;
        end

        function [T,unmatched] = table_with_labs(this, varargin)
            if ~isempty(this.table_with_labs_internal_)
                T = this.table_with_labs_internal_;
                return
            end

            T = table(this);

            % add labs vars
            date_labs = NaT(size(T, 1), 1);
            T = addvars(T, date_labs);  % for Tlabs.DateOfPETScan
            glc_mg_dL = nan(size(T, 1), 1);
            T = addvars(T, glc_mg_dL);
            PulseOx = nan(size(T, 1), 1);
            T = addvars(T, PulseOx);
            hct = nan(size(T, 1), 1);
            T = addvars(T, hct);
            hgb = nan(size(T, 1), 1);
            T = addvars(T, hgb);

            % assign labs vars
            labs = mlvg.Laboratory();
            Tlabs = table(labs);
            unmatched = 0;
            for rec_idx = 1:length(Tlabs.record_id)
                try
                    new_pos = find(T.record_id == Tlabs.record_id(rec_idx));
                    if isempty(new_pos)
                        unmatched = unmatched + 1;
                        continue
                    end
                    if ~isscalar(new_pos)
                        pet_dates = T.pet_metab_last_date(new_pos);
                        lab_date = Tlabs.DateOfPETScan(rec_idx);
                        separation_dates = abs(days(pet_dates - lab_date));
                        [~, sep_idx] = min(separation_dates);
                        new_pos = new_pos(sep_idx);
                    end
                    T.date_labs(new_pos) = Tlabs.DateOfPETScan(rec_idx);
                    T.glc_mg_dL(new_pos) = Tlabs.glc_mg_dL(rec_idx);
                    T.PulseOx(new_pos) = Tlabs.PulseOx(rec_idx);
                    T.hct(new_pos) = Tlabs.hct(rec_idx);
                    T.hgb(new_pos) = Tlabs.hgb(rec_idx);
                catch ME
                    handwarning(ME)
                end
            end
            T.date_labs.Format = "yyyy-MM-dd";

            % manage missing hct, hgb
            % T = this.impute_hgb(T);

            this.table_with_labs_internal_ = T;
        end
    end

    methods (Static)
        function PaO2 = EllisSeveringhaus(SaO2)
            %% Roger K. Ellis.  https://www.physiology.org/doi/10.1152/jappl.1989.67.2.902
            %  See also:  John W. Severinghaus. Simple, accurate equations for
            %  human blood O2 dissociation computations. J. Appl. Physiol:
            %  Respirat. Environ. Exercise Physiol. 46(3):599-602, 1979.
            %  revisions, 1999, 2002, 2007.
            %  Returns mm Hg
            %  
            %  Check match with Severinghaus' empirical with:
            %      S = 0.5:.01:0.99;  plot(S, mlvg.Clinical.EllisSeveringhaus(S))  
            
            % manage nan and fractions
            isnice = ~isnan(SaO2);
            if all(SaO2(isnice) > 1)
                SaO2 = SaO2 / 100;
            end
            assert(all(0.49 < SaO2(isnice) & SaO2(isnice) <= 1))

            % Ellis' equation
            A = 11700 ./ (1 ./ SaO2 - 1);
            B = sqrt(50^3 + A.^2);
            PaO2 = (A + B).^(1/3) + (A - B).^(1/3);
            PaO2 = real(PaO2);  % mm Hg

            % practical bound at 105 mm Hg
            if any(SaO2 > 0.98)
                high = SaO2 > 0.98;
                PaO2(high) = 105;
                return
            end
        end
        
        function CaO2 = Marino_eq_10_4(hgb, SaO2)
            %% Marino's ICU Book.  CaO2 = 1.34 [Hb] SaO2 + 0.03 PaO2 (10.4). 
            %  The O2 content in arterial blood (CaO2) is determined by combining 
            %  equations 10.2, hgb-bound oxygen, and 10.3, dissolved oxygen, and 
            %  inserting the SO2 and PO2 of arterial blood (SaO2 and PaO2).
            %  Returns mL/dL
            
            %  Hemoglobin-Bound Oxygen.
            %  [Hb] is the concentration of hemoglobin in g/dL (grams per 100 mL), 
            %  1.34 is the O2 binding capacity of hemoglobin, in mL/g 
            %  (i.e., one gram of hemoglobin will bind 1.34 mL of O2 when fully 
            %  saturated), and SO2 is the O2 saturation, expressed as a ratio.

            assert(isnumeric(hgb))
            assert(all(5 < hgb & hgb < 25))

            %  Dissolved Oxygen.
            %  Oxygen does not dissolve readily in plasma (which is why hemoglobin 
            %  is needed as a carrier molecule). The solubility of O2 in plasma is 
            %  temperature-dependent, and varies inversely with a change in body 
            %  temperature. At a normal body temperature (37°C), each increment in 
            %  PO2 of 1 mm Hg will increase the concentration of dissolved O2 by 
            %  0.03 mL/L (4). This relationship is expressed as a solubility 
            %  coefficient of 0.03 mL/L/mm Hg. 

            % manage nan and fractions
            isnice = ~isnan(SaO2);
            if all(SaO2(isnice) > 1)
                SaO2 = SaO2 / 100;
            end
            assert(all(0.49 < SaO2(isnice) & SaO2(isnice) <= 1))

            PaO2 = mlvg.Clinical.EllisSeveringhaus(SaO2);

            CaO2 = 1.34 * hgb .* SaO2 + 0.003 .* PaO2;  % mL/dL dissolved in arterial blood

            %  As shown in Table 10.1, the normal arterial O2 content is 20 mL/dL (or 200 mL/L), 
            %  and only 1.5% (0.3 mL/dL) represents dissolved O2.  Note also that the total 
            %  volume of O2 in arterial blood is less than half the volume of O2 in venous blood (!). 
            %  This is a reflection of the uneven distribution of blood volume in the circulatory 
            %  system, with 75% of the volume in the veins.
        end
    
        function mumol_mL = molar_density_glc(mg_dL)
            % Glucose (μmol/mL) = Glucose (mg/dL) × \frac{1 mmol​}{180.16 mg} × \frac{1000 μmol​}{1 mmol} × \frac{1 dL}{100 mL}​
            % Glucose (μmol/mL) = Glucose (mg/dL) × 0.0555

            mumol_mL = mg_dL * 0.0555;  % divide by 18.0
        end

        function mumol_mL = molar_density_O2(mL_dL)
            % To convert from mL O₂/dL to μmoles/mL:
            % CaO2​ (μmoles/mL) = CaO2​ (mL/dL) × \frac{1000 μmoles}{22.4 mL} × \frac{1 dL}{100 mL}
            % CaO2​ (μmoles/mL) = CaO2​ (mL/dL) × 0.446

            mumol_mL = mL_dL * 0.446;
        end
    end

    %% PRIVATE

    properties (Access = private)
        table_with_labs_internal_
        table_internal_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
