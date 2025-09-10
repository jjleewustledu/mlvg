classdef Laboratory < handle
    %% line1
    %  line2
    %  
    %  Created 06-Sep-2025 13:48:19 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.

        
    methods
        function this = Laboratory(xlsx_fqfn)
            arguments
                xlsx_fqfn {mustBeFile} = ...
                    fullfile(getenv("HOME"), "MATLAB-Drive", "mlvg", "data", ...
                    "Vision PET scans with pulse ox and glucose data.xlsx")
                % see also Kim Casey's emails from 20250905
            end
            
            T = readtable(xlsx_fqfn);            

            % simplify variables
            T.record_id = uint64(T.record_id);
            pulse_oxs = [T.PulseOx_CO, T.PulseOx_O2, T.PulseOx_O2_1, T.PulseOx_H2O, T.PulseOx_H2O_1];
            PulseOx = median(pulse_oxs, 2, "omitnan");
            PulseOx(isnan(PulseOx)) = median(PulseOx, 1, "omitnan");
            T = addvars(T, PulseOx);
            selection = ~isnan(T.bloodSugarRetakenCloserToFDGScan);
            T.iv_fdg_blood_sugar_level_wmh(selection) = T.bloodSugarRetakenCloserToFDGScan(selection);
            T = renamevars(T, "iv_fdg_blood_sugar_level_wmh", "glc_mg_dL");
            T = removevars(T, [ ...
                "PulseOx_CO", "PulseOx_O2", "PulseOx_O2_1", "PulseOx_H2O", "PulseOx_H2O_1", ...
                "ScanStartTimeCO", "ScanStartTimeO2", "ScanStartTimeO2_1", "ScanStartTimeH2O", "H2ORepeated_StartScanTime", ...
                "bloodSugarRetakenCloserToFDGScan"]);

            this.table_internal_ = T;
        end

        function T = table(this, varargin)
            T = this.table_internal_;
        end
    end

    methods (Static)
        function hgb = Nijboer_2007(hct)
            % Myth or Reality: Hematocrit and Hemoglobin Differ in Trauma
            % Nijboer, Johanna M. M. MD; van der Horst, Iwan C. C. MD, PhD; Hendriks, Herman G. D. MD, PhD; 
            % ten Duis, Hendrik-Jan MD, PhD; Nijsten, Maarten W. N. MD, PhD
            % The Journal of Trauma: Injury, Infection, and Critical Care 62(5):p 1310-1312, May 2007.
            % DOI: 10.1097/TA.0b013e3180341f54
            %
            % Arg:  hct := "m"|"f"|fraction|percent
            % Return:  hgb in g/dL

            % manage text
            if istext(hct)
                hgb = nan(size(hct));
                ismale = strcmpi(hct, "m");
                isfem = strcmpi(hct, "f");
                if any(ismale)
                    hgb(ismale) = 15.50;  % g/dL
                end
                if any(isfem)
                    hgb(isfem) = 13.75;  % g/dL
                end
                return
            else
                assert(isnumeric(hct))
            end

            % manage nan and fractions
            isnice = ~isnan(hct);
            if all(hct(isnice) < 1)
                hct = 100 * hct;
            end
            assert(all(hct(isnice) < 100))

            hgb = 0.334 * hct;  % g/dL
        end

    end

    %% PRIVATE

    properties (Access = private)
        table_internal_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
