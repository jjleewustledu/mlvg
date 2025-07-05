classdef Idif2024Par
    %% See also PycharmProjects/dynesty/idif204.
    %  1 parcel requires < 0.35 GB memory, 309 instances of multiprocessing.Pool, <1 h for 100 nlive, ~24 h for 1000 nlive.
    %  For 100 nlive:  request nodes with 11 GB memory, 32 cores, 10 h.  
    %  For 1000 nlive:  request nodes with 11 GB memory, 32 cores, 240 h.
    %  
    %  Created 29-Feb-2024 00:16:53 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 9.12.0.2327980 (R2022a) Update 7 for MACI64.  Copyright 2024 John J. Lee.
    
    methods (Static)
        function j = serial_tiny_test()
            mlvg.CHPC3.propcluster_tiny()
            
            c = parcluster;
            disp(c.AdditionalProperties)
            % Submit job to query where MATLAB is running on the cluster
            j = c.batch(@pwd, 1, {}, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false);

            % Query job for state
            pause(60)
            j.State

            % If state is finished, fetch the results
            if contains(j.State, 'finished')
                j.fetchOutputs{:}
            end
        end
        function j = serial_main4_test()
            mlvg.CHPC3.propcluster()

            c = parcluster;
            disp(c.AdditionalProperties)
            j = c.batch(@mlvg.Idif2024Par.main, 1, {'version', 4}, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false);
        end

        function main(opts)
            arguments
                opts.machine {mustBeTextScalar} = "cluster"
                opts.version {mustBeNumeric} = 4
            end
            opts.version = string(opts.version);
            opts.version = strrep(opts.version, ".", "_");

            switch convertStringsToChars(opts.machine)
                case {'login3', 'chpc', 'cluster'}
                    dyn_bin = "/home/jjlee/miniconda/envs/dynesty/bin";
                case {'vglab2', 'linux1', 'linux2', 'pascal'}
                    dyn_bin = "/data/nil-bluearc/raichle/jjlee/anaconda3/envs/dynesty/bin";
                otherwise
                    error("mlvg:ValueError", stackstr())
            end
            setenv("PATH", dyn_bin+":"+getenv("PATH"));

            main_py = fullfile( ...
                getenv("HOME"), "PycharmProjects", "dynesty", "idif2024", "main"+opts.version+".py");
            pyrunfile(main_py)
        end

        
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
