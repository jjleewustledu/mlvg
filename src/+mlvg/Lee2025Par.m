classdef Lee2025Par < handle & mlvg.Lee2025
    %% Use cluster_call:  ho steps 1,4 first; then others steps 1, 4; then all step 5.
    %  
    %  Created 04-Jul-2025 01:08:45 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.
    

    methods
        function this = Lee2025Par(varargin)
            this = this@mlvg.Lee2025(varargin{:});
        end
    end

    methods (Static)
        function [j,c,msg,id] = cluster_call(globbing_mat, opts)
            %% for clusters running Matlab parallel server

            arguments
                globbing_mat {mustBeFile} = ...
                    fullfile( ...
                    getenv("SINGULARITY_HOME"), "CCIR_01211", "mlvg_Lee2025Par_globbing_oo.mat")
                opts.globbing_var = "globbed"
                opts.selection_indices double = []  % total ~ 1:58 for ho, 1:69 for co, 1:112 for oo
                opts.Ncol {mustBeInteger} = 6
                opts.method {mustBeTextScalar} = "do_make_input_func"
                opts.reference_tracer {mustBeTextScalar} = "ho"
                opts.steps {mustBeNumericOrLogical} = 4
                opts.account {mustBeTextScalar} = "manu_goyal"
            end
            ld = load(globbing_mat);
            globbed = convertCharsToStrings(ld.(opts.globbing_var));
            globbed = asrow(globbed);
            if ~isempty(opts.selection_indices)
                globbed = globbed(opts.selection_indices);
            end

            % pad and reshape globbed
            Nrow = ceil(numel(globbed)/opts.Ncol);
            padding = repmat("", [1, opts.Ncol*Nrow - numel(globbed)]);
            globbed = [globbed, padding];
            globbed = reshape(globbed, Nrow, opts.Ncol);
            fprintf("%s:globbed:\n", stackstr())
            disp(size(globbed))
            disp(ascol(globbed))

            % contact cluster slurm

            warning('off', 'MATLAB:legacy:batchSyntax');
            warning('off', 'parallel:convenience:BatchFunctionNestedCellArray');
            warning('off', 'MATLAB:TooManyInputs');

            c = mlvg.CHPC3.propcluster(opts.account, mempercpu='32gb', walltime='8:00:00');
            disp(c.AdditionalProperties)
            for irow = 1:Nrow
                try
                    j = c.batch( ...
                        @mlvg.Lee2025Par.par_call_ifk, ...
                        1, ...
                        {globbed(irow, :), ...
                        'method', opts.method, 'steps', opts.steps, 'reference_tracer', opts.reference_tracer}, ...
                        'Pool', opts.Ncol, ...
                        'CurrentFolder', '/scratch/jjlee/Singularity/CCIR_01211', ...
                        'AutoAddClientPath', false);
                catch ME
                    handwarning(ME)
                end
            end

            [msg,id] = lastwarn();
        end

        function durations = par_call_ifk(nii, opts)
            arguments
                nii {mustBeText} = "sourcedata/sub-108007/ses-20210219145054/pet/sub-108007_ses-20210219145054_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz"
                opts.out_dir {mustBeFolder} = "/scratch/jjlee/Singularity/CCIR_01211"
                opts.method {mustBeTextScalar} = "do_make_input_func"
                opts.reference_tracer {mustBeTextScalar} = "ho"
                opts.steps {mustBeNumericOrLogical} = 1
            end
            
            nii = nii(~arrayfun(@isempty, nii));  % Remove empty cells
            durations = nan(1, length(nii));

            parfor sidx = 1:length(nii)

                tic;
            
                % setup
                mlvg.CHPC3.setenvs();

                try
                    % construct & call
                    lp = mlvg.Lee2025Par(nii(sidx), out_dir=opts.out_dir); %#ok<PFBNS>
                    call_ifk(lp, method=opts.method, steps=opts.steps, reference_tracer=opts.reference_tracer);
                catch ME
                    handwarning(ME)
                end

                durations(sidx) = toc;
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
