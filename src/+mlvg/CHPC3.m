classdef CHPC3
    %% Provides support functions for using the Matlab Parallel Server at CHPC3.
    %  
    %  Created 07-Apr-2022 16:13:53 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.10.0.1851785 (R2021a) Update 6 for MACI64.  Copyright 2022 John J. Lee.
    
    methods (Static)
        function clean_tempdir()
            try
                deleteExisting(fullfile(tempdir, '*.nii*'));
                deleteExisting(fullfile(tempdir, '*.mat'));
                deleteExisting(fullfile(tempdir, '*.save'));
                deleteExisting(fullfile(tempdir, 'tp*'));
            catch ME
                disp(ME)
            end
        end

        function c = propcluster(account_name, opts)
            arguments
                account_name = 'manu_goyal'  % 'aristeidis_sotiras' 'joshua_shimony' 'manu_goyal' 'john_lee'
                opts.partition = 'tier1_cpu'  % 'tier2_cpu' 'tier1_cpu'
                opts.mempercpu {mustBeTextScalar} = '64gb'
                opts.walltime {mustBeTextScalar} = '01:00:00'
            end
            if strcmp(account_name, 'aristeidis_sotiras')
                opts.partition = 'tier2_cpu';
            end
            account_name = convertStringsToChars(account_name);
            opts.partition = convertStringsToChars(opts.partition);

            c = parcluster;
            c.AdditionalProperties.AccountName = account_name;
            c.AdditionalProperties.AdditionalSubmitArgs = sprintf('--account=%s', account_name);
            c.AdditionalProperties.ClusterHost = 'login3.chpc.wustl.edu';
            c.AdditionalProperties.EmailAddress = '';
            c.AdditionalProperties.EnableDebug = false;
            c.AdditionalProperties.GpusPerNode = 0;
            c.AdditionalProperties.MemPerCPU = opts.mempercpu;
            % c.AdditionalProperties.Node = '';
            c.AdditionalProperties.Partition = opts.partition;
            c.AdditionalProperties.RemoteJobStorageLocation = '/home/jjlee/.matlab/3p_cluster_jobs/chpc/twistor.attlocal.net.dhcp.wustl.edu/R2024b/nonshared';
            c.AdditionalProperties.UseIdentityFile = false;
            c.AdditionalProperties.UseSmpd = false;
            c.AdditionalProperties.Username = 'jjlee';
            c.AdditionalProperties.WallTime = opts.walltime;
            disp(c.AdditionalProperties)

            warning('off', 'MATLAB:legacy:batchSyntax');
            warning('off', 'parallel:convenience:BatchFunctionNestedCellArray');
            warning('off', 'MATLAB:TooManyInputs');
            warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary')
        end

        function propcluster_tiny()
            c = parcluster;
            c.AdditionalProperties.EmailAddress = '';
            c.AdditionalProperties.EnableDebug = 0;
            c.AdditionalProperties.GpusPerNode = 0;
            c.AdditionalProperties.MemUsage = '32000'; % in MB
            c.AdditionalProperties.Node = 1;
            c.AdditionalProperties.Partition = 'tier2_cpu';
            c.AdditionalProperties.AdditionalSubmitArgs = '--account=aristeidis_sotiras';
            c.AdditionalProperties.WallTime = '01:00:00';
            c.saveProfile
            disp(c.AdditionalProperties)
        end
        function propcluster_free()
            c = parcluster;
            c.AdditionalProperties.EmailAddress = '';
            c.AdditionalProperties.EnableDebug = 0;
            c.AdditionalProperties.GpusPerNode = 0;
            c.AdditionalProperties.MemUsage = '4000'; % in MB
            c.AdditionalProperties.Node = 1;
            c.AdditionalProperties.Partition = 'free';
            c.AdditionalProperties.AdditionalSubmitArgs = '';
            c.AdditionalProperties.WallTime = '01:00:00';
            c.saveProfile
            disp(c.AdditionalProperties)
        end
        function setenvs()
            if ~isInParallelWorker()
                return
            end

            setenv('DEBUG_SETENVS', '1')
            setenv('FSLDIR', '/scratch/jjlee/fsl')
            setenv('TMPDIR', '/scratch/jjlee/tmp')
            setenv('APPTAINER_HOME', '/scratch/jjlee/Singularity')
            setenv('SINGULARITY_HOME', '/scratch/jjlee/Singularity')
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
