classdef RedCap
    %% line1
    %  line2
    %  
    %  Created 06-Apr-2023 18:20:32 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.14.0.2206163 (R2023a) for MACI64.  Copyright 2023 John J. Lee.
    
    methods (Static)
        function build_table()
            load('O15FDGVenous.mat');
            T0 = O15FDGVenous;
            for c = 1:260
                col = T0{:,c};
                if isstring(col)
                    if all(col == ""); fprintf('empty string col -> %i', c); end
                end
                if isnumeric(col)
                    if all(isnan(col)); fprintf('nan col -> %i', c); end
                end
                if isdatetime(col)
                    if all(isnat(col)); fprintf('nat col -> %i', c); end
                end
            end

            T = table( ...
                T0.RecordID, T0.EventName, ...
                T0.PleaseIncludeAnyNotesAboutO2CountsHere, T0.DateO2CountsWereAcquired, T0.Complete, ...
                T0.PleaseIncludeAnyNotesAboutFDGCountsHere, T0.DateFDGCountsWereAcquired, T0.Complete1 ...
                );
            T.Properties.VariableNames = { ...
                'RecordID', 'EventName', ...
                'O15Notes', 'O15Date', 'O15Complete', ...
                'FDGNotes', 'FDGDate', 'FDGComplete'};
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
