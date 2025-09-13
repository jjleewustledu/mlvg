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
        function [j,c,msg,id] = cluster_test(globbing_mat, opts)
            %% for clusters running Matlab parallel server
            %
            %  E.g.:
            %  [j,c,msg,id] = mlvg.Lee2025Par.cluster_test()
            %  j_output = fetchOutputs(j)
            %  fprintf("%s\n", j_output{1})
            %
            % fprintf("%s\n", j_output{1}) % ================== without mlenv.sh ==================
            % ARCH=glnxa64
            % AUTOMOUNT_MAP=
            % BASEMATLABPATH=
            % DEBUG_SETENVS=1
            % DISPLAY=
            % GFORTRAN_UNBUFFERED_PRECONNECTED=y
            % HOME=/home/jjlee
            % HYDRA_USER_PROVIDED_BINDING=1
            % ICU_TIMEZONE_FILES_DIR=/matlab-local/R2024b/bin/icutzdata
            % KMP_BLOCKTIME=1
            % KMP_HANDLE_SIGNALS=0
            % KMP_INIT_AT_FORK=false
            % KMP_STACKSIZE=512k
            % LANG=C
            % LD_LIBRARY_PATH=/matlab-local/R2024b/sys/opengl/lib/glnxa64:/matlab-local/R2024b/sys/os/glnxa64:/matlab-local/R2024b/bin/glnxa64:/home/jjlee/.MathWorks/t_C_7MXF5pBGVIY_VqFQJOHvHuJz_gWW_muyxourG40/bin/glnxa64:/matlab-local/R2024b/extern/lib/glnxa64:/matlab-local/R2024b/cefclient/sys/os/glnxa64:/matlab-local/R2024b/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/matlab-local/R2024b/sys/java/jre/glnxa64/jre/lib/amd64/server
            % LD_PRELOAD=
            % MATLAB_MEM_MGR_VALUE=
            % MDCE_FORCE_MPI_OPTION=mpich
            % MDCE_PARALLEL=1
            % MEMKIND_HEAP_MANAGER=TBB
            % MKL_DOMAIN_NUM_THREADS=
            % MKL_NUM_THREADS=
            % MLM_WEB_LICENSE=false
            % MPIR_CVAR_CH3_INTERFACE_HOSTNAME=node16
            % MPI_LOCALNRANKS=2
            % MPI_LOCALRANKID=0
            % OLDPWD=/matlab-local/R2024b
            % OMPI_MCA_orte_precondition_transports=0058c18700000000-0058c18700000000
            % OSG_LD_LIBRARY_PATH=/matlab-local/R2024b/sys/openscenegraph/lib/glnxa64
            % PARALLEL_SERVER_CMR=/matlab-local/R2024b
            % PARALLEL_SERVER_DEBUG=false
            % PARALLEL_SERVER_DECODE_FUNCTION=parallel.cluster.generic.communicatingDecodeFcn
            % PARALLEL_SERVER_JOB_LOCATION=Job7
            % PARALLEL_SERVER_MATLAB_ARGS=-parallel
            % PARALLEL_SERVER_MATLAB_EXE=/matlab-local/R2024b/bin/worker
            % PARALLEL_SERVER_NUM_THREADS=1
            % PARALLEL_SERVER_STORAGE_CONSTRUCTOR=makeFileStorageObject
            % PARALLEL_SERVER_STORAGE_LOCATION=/home/jjlee/.matlab/3p_cluster_jobs/chpc/twistor.attlocal.net.dhcp.wustl.edu/R2024b/nonshared
            % PARALLEL_SERVER_TOTAL_TASKS=2
            % PMIX_BFROP_BUFFER_TYPE=PMIX_BFROP_BUFFER_NON_DESC
            % PMIX_DSTORE_21_BASE_PATH=/var/spool/slurmd/pmix.5816711.0//pmix_dstor_ds21_1050522
            % PMIX_DSTORE_ESH_BASE_PATH=/var/spool/slurmd/pmix.5816711.0//pmix_dstor_ds12_1050522
            % PMIX_GDS_MODULE=ds21,ds12,hash
            % PMIX_HOSTNAME=node16
            % PMIX_NAMESPACE=slurm.pmix.5816711.0
            % PMIX_RANK=0
            % PMIX_SECURITY_MODE=munge,native
            % PMIX_SERVER_TMPDIR=/var/spool/slurmd/pmix.5816711.0/
            % PMIX_SERVER_URI21=pmix-server.1050522;tcp4://127.0.0.1:60043
            % PMIX_SERVER_URI2=pmix-server.1050522;tcp4://127.0.0.1:60043
            % PMIX_SERVER_URI3=pmix-server.1050522;tcp4://127.0.0.1:60043
            % PMIX_SERVER_URI41=pmix-server.1050522;tcp4://127.0.0.1:60043
            % PMIX_SERVER_URI4=pmix-server.1050522;tcp4://127.0.0.1:60043
            % PMIX_SYSTEM_TMPDIR=/tmp
            % PMIX_VERSION=4.2.9
            % PMI_FD=6
            % PMI_RANK=0
            % PMI_SIZE=2
            % PRE_LD_PRELOAD=
            % PWD=/scratch/jjlee/Singularity/CCIR_01211
            % SHELL=/bin/bash
            % SHLVL=1
            % SLURMD_DEBUG=2
            % SLURMD_NODENAME=node16
            % SLURM_CLUSTER_NAME=chpc3
            % SLURM_CONF=/etc/slurm/slurm.conf
            % SLURM_CPUS_ON_NODE=2
            % SLURM_CPUS_PER_TASK=1
            % SLURM_CPU_BIND=quiet,mask_cpu:0x30000000
            % SLURM_CPU_BIND_LIST=0x30000000
            % SLURM_CPU_BIND_TYPE=mask_cpu:
            % SLURM_CPU_BIND_VERBOSE=quiet
            % SLURM_DISTRIBUTION=cyclic
            % SLURM_EXPORT_ENV=NONE
            % SLURM_GET_USER_ENV=1
            % SLURM_GTIDS=0
            % SLURM_JOBID=5816711
            % SLURM_JOB_ACCOUNT=manu_goyal
            % SLURM_JOB_CPUS_PER_NODE=2
            % SLURM_JOB_END_TIME=1753988347
            % SLURM_JOB_GID=1000070
            % SLURM_JOB_ID=5816711
            % SLURM_JOB_NAME=MATLAB_R2024b_Job7
            % SLURM_JOB_NODELIST=node16
            % SLURM_JOB_NUM_NODES=1
            % SLURM_JOB_PARTITION=tier1_cpu
            % SLURM_JOB_QOS=normal
            % SLURM_JOB_START_TIME=1753987747
            % SLURM_JOB_UID=1002720
            % SLURM_JOB_USER=jjlee
            % SLURM_LAUNCH_NODE_IPADDR=10.1.1.16
            % SLURM_LOCALID=0
            % SLURM_MEM_BIND=quiet,none
            % SLURM_MEM_BIND_LIST=
            % SLURM_MEM_BIND_TYPE=none
            % SLURM_MEM_BIND_VERBOSE=quiet
            % SLURM_MEM_PER_CPU=16384
            % SLURM_MPI_TYPE=pmix_v4
            % SLURM_NNODES=1
            % SLURM_NODEID=0
            % SLURM_NODELIST=node16
            % SLURM_NPROCS=1
            % SLURM_NTASKS=1
            % SLURM_NTASKS_PER_CORE=1
            % SLURM_OOM_KILL_STEP=0
            % SLURM_PMIXP_ABORT_AGENT_PORT=41297
            % SLURM_PMIX_MAPPING_SERV=(vector,(0,1,1))
            % SLURM_PRIO_PROCESS=0
            % SLURM_PROCID=0
            % SLURM_SRUN_COMM_HOST=10.1.1.16
            % SLURM_SRUN_COMM_PORT=42447
            % SLURM_STEPID=0
            % SLURM_STEP_ID=0
            % SLURM_STEP_LAUNCHER_PORT=42447
            % SLURM_STEP_NODELIST=node16
            % SLURM_STEP_NUM_NODES=1
            % SLURM_STEP_NUM_TASKS=1
            % SLURM_STEP_TASKS_PER_NODE=1
            % SLURM_SUBMIT_DIR=/home/jjlee
            % SLURM_SUBMIT_HOST=login01.cluster
            % SLURM_TASKS_PER_NODE=1
            % SLURM_TASK_PID=1050617
            % SLURM_TOPOLOGY_ADDR=node16
            % SLURM_TOPOLOGY_ADDR_PATTERN=node
            % SLURM_TRES_PER_TASK=cpu=1
            % SLURM_UMASK=0022
            % TMPDIR=/tmp
            % TOOLBOX=/matlab-local/R2024b/toolbox
            % TZ=America/Chicago
            % USER=jjlee
            % XFILESEARCHPATH=/matlab-local/R2024b/sys/java/jre/glnxa64/jre/lib/locale/%L/%T/%N%S:
            % ZES_ENABLE_SYSMAN=1
            % _=/usr/bin/env
            %
            % fprintf("%s\n", j_output{1})  % ================== with mlenv.sh =======================
            %  *v*)
            %  *v*x*)
            %  *x*)
            %  ;;
            %  ;;
            %  ;;
            %  __lmod_my_status=$?;
            %  __lmod_sh_dbg='v'
            %  __lmod_sh_dbg='vx'
            %  __lmod_sh_dbg='x'
            %  case "$-" in
            %  echo "Shell debugging restarted" 1>&2;
            %  echo "Shell debugging temporarily silenced: export LMOD_SH_DBG_ON=1 for Lmod's output" 1>&2;
            %  esac;
            %  eval "$($LMOD_CMD shell "$@")" && eval "$(${LMOD_SETTARG_CMD:-:} -s sh)";
            %  eval ${which_declare} ) | /usr/bin/which --tty-only --read-alias --read-functions --show-tilde --show-dot $@
            %  fi;
            %  fi;
            %  fi;
            %  if [ -n "${__lmod_sh_dbg:-}" ]; then
            %  if [ -n "${__lmod_sh_dbg:-}" ]; then
            %  return $__lmod_my_status
            %  set +$__lmod_sh_dbg;
            %  set -$__lmod_sh_dbg;
            %  unset __lmod_sh_dbg;
            % ADNI_FDG=/scratch/jjlee/Singularity/ADNI
            % ADNI_HOME=/scratch/jjlee/Singularity/ADNI
            % ARCH=glnxa64
            % AUTOMOUNT_MAP=
            % BASEMATLABPATH=
            % BASH_ENV=/usr/share/lmod/lmod/init/bash
            % BASH_FUNC_ml%%=() {  eval "$($LMOD_DIR/ml_cmd "$@")"
            % BASH_FUNC_module%%=() {  if [ -z "${LMOD_SH_DBG_ON+x}" ]; then
            % BASH_FUNC_which%%=() {  ( alias;
            % CCIR_RAD_MEASUREMENTS_DIR=/scratch/jjlee/Singularity/CCIR_RAD_MEASUREMENTS
            % CONDA=/home/jjlee/miniconda
            % CONTAINER_PATH=/containers/fsl-6.0.7.18-amd64.sif
            % CONTAINER_PATH_ADDNOISETOIMAGE=/containers/rocky8_mod.sif
            % CONTAINER_PATH_AFNI=/containers/afni-25.1.03.sif
            % CONTAINER_PATH_ANTS=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSAFFINEINITIALIZER=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSAI=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSALIGNORIGIN=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSAPPLYTRANSFORMS=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSAPPLYTRANSFORMSTOPOINTS=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSASLPROCESSINGSH=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSATROPOSN4SH=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSBOLDNETWORKANALYSISR=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSBRAINEXTRACTIONSH=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSCORTICALTHICKNESSSH=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSINTEGRATEVECTORFIELD=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSINTEGRATEVELOCITYFIELD=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSINTERMODALITYINTRASUBJECTSH=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSINTRODUCTIONSH=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSJACOBIAN=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSJOINTFUSION=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSJOINTLABELFUSIONSH=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSJOINTTENSORFUSION=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSLANDMARKBASEDTRANSFORMINITIALIZER=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSLAPLACIANBOUNDARYCONDITIONR=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSLONGITUDINALCORTICALTHICKNESSSH=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSMOTIONCORR=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSMOTIONCORRDIFFUSIONDIRECTION=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSMOTIONCORRSTATS=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSMULTIVARIATETEMPLATECONSTRUCTION2SH=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSMULTIVARIATETEMPLATECONSTRUCTIONSH=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSNETWORKANALYSISR=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSNEUROIMAGINGBATTERY=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSPEXECSH=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSREGISTRATION=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSREGISTRATIONSYNQUICKSH=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSREGISTRATIONSYNSH=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSSLICEREGULARIZEDREGISTRATION=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSTRANSFORMINFO=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSUSEDEFORMATIONFIELDTOGETAFFINETRANSFORM=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSUSELANDMARKIMAGESTOGETAFFINETRANSFORM=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSUSELANDMARKIMAGESTOGETBSPLINEDISPLACEMENTFIELD=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ANTSUTILITIESTESTING=/containers/rocky8_mod.sif
            % CONTAINER_PATH_ATROPOS=/containers/rocky8_mod.sif
            % CONTAINER_PATH_AVERAGEAFFINETRANSFORM=/containers/rocky8_mod.sif
            % CONTAINER_PATH_AVERAGEAFFINETRANSFORMNORIGID=/containers/rocky8_mod.sif
            % CONTAINER_PATH_AVERAGEIMAGES=/containers/rocky8_mod.sif
            % CONTAINER_PATH_AVERAGETENSORIMAGES=/containers/rocky8_mod.sif
            % CONTAINER_PATH_BBREGISTER=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_CLUSTERIMAGESTATISTICS=/containers/rocky8_mod.sif
            % CONTAINER_PATH_COMPARETWOTRANSFORMS=/containers/rocky8_mod.sif
            % CONTAINER_PATH_COMPOSEMULTITRANSFORM=/containers/rocky8_mod.sif
            % CONTAINER_PATH_COMPOSITETRANSFORMUTIL=/containers/rocky8_mod.sif
            % CONTAINER_PATH_CONVERTIMAGE=/containers/rocky8_mod.sif
            % CONTAINER_PATH_CONVERTIMAGEPIXELTYPE=/containers/rocky8_mod.sif
            % CONTAINER_PATH_CONVERTINPUTIMAGEPIXELTYPETOFLOAT=/containers/rocky8_mod.sif
            % CONTAINER_PATH_CONVERTSCALARIMAGETORGB=/containers/rocky8_mod.sif
            % CONTAINER_PATH_CONVERTTOJPG=/containers/rocky8_mod.sif
            % CONTAINER_PATH_CONVERTTRANSFORMFILE=/containers/rocky8_mod.sif
            % CONTAINER_PATH_COPYIMAGEHEADERINFORMATION=/containers/rocky8_mod.sif
            % CONTAINER_PATH_CPP=/containers/rocky8_mod.sif
            % CONTAINER_PATH_CREATEDISPLACEMENTFIELD=/containers/rocky8_mod.sif
            % CONTAINER_PATH_CREATEDTICOHORT=/containers/rocky8_mod.sif
            % CONTAINER_PATH_CREATEIMAGE=/containers/rocky8_mod.sif
            % CONTAINER_PATH_CREATEJACOBIANDETERMINANTIMAGE=/containers/rocky8_mod.sif
            % CONTAINER_PATH_CREATETILEDMOSAIC=/containers/rocky8_mod.sif
            % CONTAINER_PATH_CREATEWARPEDGRIDIMAGE=/containers/rocky8_mod.sif
            % CONTAINER_PATH_C__=/containers/rocky8_mod.sif
            % CONTAINER_PATH_DENOISEIMAGE=/containers/rocky8_mod.sif
            % CONTAINER_PATH_DENRRD=/containers/rocky8_mod.sif
            % CONTAINER_PATH_EXTRACTREGIONFROMIMAGE=/containers/rocky8_mod.sif
            % CONTAINER_PATH_EXTRACTREGIONFROMIMAGEBYMASK=/containers/rocky8_mod.sif
            % CONTAINER_PATH_EXTRACTSLICEFROMIMAGE=/containers/rocky8_mod.sif
            % CONTAINER_PATH_FITBSPLINETOPOINTS=/containers/rocky8_mod.sif
            % CONTAINER_PATH_FOOSH=/containers/rocky8_mod.sif
            % CONTAINER_PATH_FREEVIEW=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_GCC=/containers/rocky8_mod.sif
            % CONTAINER_PATH_GCC_AR=/containers/rocky8_mod.sif
            % CONTAINER_PATH_GCC_NM=/containers/rocky8_mod.sif
            % CONTAINER_PATH_GCC_RANLIB=/containers/rocky8_mod.sif
            % CONTAINER_PATH_GCOV=/containers/rocky8_mod.sif
            % CONTAINER_PATH_GCOV_DUMP=/containers/rocky8_mod.sif
            % CONTAINER_PATH_GCOV_TOOL=/containers/rocky8_mod.sif
            % CONTAINER_PATH_GETCONNECTEDCOMPONENTSFEATUREIMAGES=/containers/rocky8_mod.sif
            % CONTAINER_PATH_GFORTRAN=/containers/rocky8_mod.sif
            % CONTAINER_PATH_G__=/containers/rocky8_mod.sif
            % CONTAINER_PATH_IMAGECOMPARE=/containers/rocky8_mod.sif
            % CONTAINER_PATH_IMAGEINTENSITYSTATISTICS=/containers/rocky8_mod.sif
            % CONTAINER_PATH_IMAGEMATH=/containers/rocky8_mod.sif
            % CONTAINER_PATH_IMAGESETSTATISTICS=/containers/rocky8_mod.sif
            % CONTAINER_PATH_IMATH=/containers/rocky8_mod.sif
            % CONTAINER_PATH_KELLYKAPOWSKI=/containers/rocky8_mod.sif
            % CONTAINER_PATH_KELLYSLATER=/containers/rocky8_mod.sif
            % CONTAINER_PATH_LABELCLUSTERSUNIQUELY=/containers/rocky8_mod.sif
            % CONTAINER_PATH_LABELGEOMETRYMEASURES=/containers/rocky8_mod.sif
            % CONTAINER_PATH_LABELOVERLAPMEASURES=/containers/rocky8_mod.sif
            % CONTAINER_PATH_LAPLACIANTHICKNESS=/containers/rocky8_mod.sif
            % CONTAINER_PATH_LESIONFILLING=/containers/rocky8_mod.sif
            % CONTAINER_PATH_LTO_DUMP=/containers/rocky8_mod.sif
            % CONTAINER_PATH_MEASUREIMAGESIMILARITY=/containers/rocky8_mod.sif
            % CONTAINER_PATH_MEASUREMINMAXMEAN=/containers/rocky8_mod.sif
            % CONTAINER_PATH_MEMORYTEST=/containers/rocky8_mod.sif
            % CONTAINER_PATH_MESAGL_WB_COMMAND=/containers/workbench-2.0.1.sif
            % CONTAINER_PATH_MESAGL_WB_VIEW=/containers/workbench-2.0.1.sif
            % CONTAINER_PATH_MRIS_ANATOMICAL_STATS=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_MRIS_CONVERT=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_MRI_APARC2ASEG=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_MRI_BINARIZE=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_MRI_CONCAT=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_MRI_CONVERT=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_MRI_COREG=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_MRI_GLMFIT=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_MRI_GLMFIT_SIM=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_MRI_INFO=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_MRI_LABEL2LABEL=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_MRI_LABEL2VOL=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_MRI_LABEL_STATS=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_MRI_MASK=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_MRI_SEGSTATS=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_MRI_SURF2VOL=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_MRI_VOL2SURF=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_MRI_VOL2VOL=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_MULTIPLYIMAGES=/containers/rocky8_mod.sif
            % CONTAINER_PATH_N3BIASFIELDCORRECTION=/containers/rocky8_mod.sif
            % CONTAINER_PATH_N4BIASFIELDCORRECTION=/containers/rocky8_mod.sif
            % CONTAINER_PATH_NONLOCALSUPERRESOLUTION=/containers/rocky8_mod.sif
            % CONTAINER_PATH_PASTEIMAGEINTOIMAGE=/containers/rocky8_mod.sif
            % CONTAINER_PATH_PERMUTEFLIPIMAGEORIENTATIONAXES=/containers/rocky8_mod.sif
            % CONTAINER_PATH_PRINTHEADER=/containers/rocky8_mod.sif
            % CONTAINER_PATH_REBASETENSORIMAGE=/containers/rocky8_mod.sif
            % CONTAINER_PATH_RECON_ALL=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_REORIENTTENSORIMAGE=/containers/rocky8_mod.sif
            % CONTAINER_PATH_RESAMPLEIMAGE=/containers/rocky8_mod.sif
            % CONTAINER_PATH_RESAMPLEIMAGEBYSPACING=/containers/rocky8_mod.sif
            % CONTAINER_PATH_RESETDIRECTION=/containers/rocky8_mod.sif
            % CONTAINER_PATH_RUN_SPM12SH=/containers/rocky8_mod.sif
            % CONTAINER_PATH_RUN_SPM12SHORIG=/containers/rocky8_mod.sif
            % CONTAINER_PATH_SCCAN=/containers/rocky8_mod.sif
            % CONTAINER_PATH_SETDIRECTIONBYMATRIX=/containers/rocky8_mod.sif
            % CONTAINER_PATH_SETORIGIN=/containers/rocky8_mod.sif
            % CONTAINER_PATH_SETSPACING=/containers/rocky8_mod.sif
            % CONTAINER_PATH_SHELL=/containers/workbench-2.0.1.sif
            % CONTAINER_PATH_SIMPLESYNREGISTRATION=/containers/rocky8_mod.sif
            % CONTAINER_PATH_SIMULATEDISPLACEMENTFIELD=/containers/rocky8_mod.sif
            % CONTAINER_PATH_SMOOTHDISPLACEMENTFIELD=/containers/rocky8_mod.sif
            % CONTAINER_PATH_SMOOTHIMAGE=/containers/rocky8_mod.sif
            % CONTAINER_PATH_SPM12_GLNXA64=/containers/rocky8_mod.sif
            % CONTAINER_PATH_SPM12_MACI64=/containers/rocky8_mod.sif
            % CONTAINER_PATH_SPM12_WIN64EXE=/containers/rocky8_mod.sif
            % CONTAINER_PATH_SPM=/containers/rocky8_mod.sif
            % CONTAINER_PATH_STACKSLICES=/containers/rocky8_mod.sif
            % CONTAINER_PATH_SUPERRESOLUTION=/containers/rocky8_mod.sif
            % CONTAINER_PATH_SURFACEBASEDSMOOTHING=/containers/rocky8_mod.sif
            % CONTAINER_PATH_SURFACECURVATURE=/containers/rocky8_mod.sif
            % CONTAINER_PATH_TEXTURECOOCCURRENCEFEATURES=/containers/rocky8_mod.sif
            % CONTAINER_PATH_TEXTURERUNLENGTHFEATURES=/containers/rocky8_mod.sif
            % CONTAINER_PATH_THRESHOLDIMAGE=/containers/rocky8_mod.sif
            % CONTAINER_PATH_TILEIMAGES=/containers/rocky8_mod.sif
            % CONTAINER_PATH_TIMESCCAN=/containers/rocky8_mod.sif
            % CONTAINER_PATH_TKREGISTER2=/containers/freesurfer-7.4.1-amd64.sif
            % CONTAINER_PATH_WAITFORPBSQJOBSPL=/containers/rocky8_mod.sif
            % CONTAINER_PATH_WAITFORSGEQJOBSPL=/containers/rocky8_mod.sif
            % CONTAINER_PATH_WAITFORSLURMJOBSPL=/containers/rocky8_mod.sif
            % CONTAINER_PATH_WAITFORXGRIDJOBSPL=/containers/rocky8_mod.sif
            % CONTAINER_PATH_WARPIMAGEMULTITRANSFORM=/containers/rocky8_mod.sif
            % CONTAINER_PATH_WARPTENSORIMAGEMULTITRANSFORM=/containers/rocky8_mod.sif
            % CONTAINER_PATH_WARPTIMESERIESIMAGEMULTITRANSFORM=/containers/rocky8_mod.sif
            % CONTAINER_PATH_WB_COMMAND=/containers/workbench-2.0.1.sif
            % CONTAINER_PATH_WB_IMPORT=/containers/workbench-2.0.1.sif
            % CONTAINER_PATH_WB_SHORTCUTS=/containers/workbench-2.0.1.sif
            % CONTAINER_PATH_WB_VIEW=/containers/workbench-2.0.1.sif
            % CONTAINER_PATH_X86_64_PC_LINUX_GNU_C__=/containers/rocky8_mod.sif
            % CONTAINER_PATH_X86_64_PC_LINUX_GNU_GCC=/containers/rocky8_mod.sif
            % CONTAINER_PATH_X86_64_PC_LINUX_GNU_GCC_1020=/containers/rocky8_mod.sif
            % CONTAINER_PATH_X86_64_PC_LINUX_GNU_GCC_AR=/containers/rocky8_mod.sif
            % CONTAINER_PATH_X86_64_PC_LINUX_GNU_GCC_NM=/containers/rocky8_mod.sif
            % CONTAINER_PATH_X86_64_PC_LINUX_GNU_GCC_RANLIB=/containers/rocky8_mod.sif
            % CONTAINER_PATH_X86_64_PC_LINUX_GNU_GFORTRAN=/containers/rocky8_mod.sif
            % CONTAINER_PATH_X86_64_PC_LINUX_GNU_G__=/containers/rocky8_mod.sif
            % DEBUG_SETENVS=1
            % DISPLAY=
            % FPATH=/usr/share/lmod/lmod/init/ksh_funcs
            % GFORTRAN_UNBUFFERED_PRECONNECTED=y
            % HCPPIPEDIR=/home/jjlee/.local/HCPpipelines-4.7.0
            % HOME=/home/jjlee
            % HYDRA_USER_PROVIDED_BINDING=1
            % ICU_TIMEZONE_FILES_DIR=/matlab-local/R2024b/bin/icutzdata
            % JULIA_DEPOT_PATH=/home/jjlee/.julia:/scratch/jjlee/shared/julia_depot
            % KMP_BLOCKTIME=1
            % KMP_HANDLE_SIGNALS=0
            % KMP_INIT_AT_FORK=false
            % KMP_STACKSIZE=512k
            % LANG=C
            % LD_LIBRARY_PATH=/usr/lib64:/matlab-local/R2024b/sys/opengl/lib/glnxa64:/matlab-local/R2024b/sys/os/glnxa64:/matlab-local/R2024b/bin/glnxa64:/home/jjlee/.MathWorks/t_C_7MXF5pBGVIY_VqFQJOHvHuJz_gWW_muyxourG40/bin/glnxa64:/matlab-local/R2024b/extern/lib/glnxa64:/matlab-local/R2024b/cefclient/sys/os/glnxa64:/matlab-local/R2024b/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/matlab-local/R2024b/sys/java/jre/glnxa64/jre/lib/amd64/server
            % LD_PRELOAD=
            % LESSOPEN=||/usr/bin/lesspipe.sh %s
            % LMOD_CMD=/usr/share/lmod/lmod/libexec/lmod
            % LMOD_DIR=/usr/share/lmod/lmod/libexec
            % LMOD_FAMILY_AFNI=afni
            % LMOD_FAMILY_AFNI_VERSION=25.1.03
            % LMOD_FAMILY_ANTS=ants
            % LMOD_FAMILY_ANTS_VERSION=2.6.0
            % LMOD_FAMILY_FREESURFER=freesurfer
            % LMOD_FAMILY_FREESURFER_VERSION=7.4.1
            % LMOD_FAMILY_FSL=fsl
            % LMOD_FAMILY_FSL_VERSION=6.0.7.18
            % LMOD_FAMILY_GCC=gcc
            % LMOD_FAMILY_GCC_VERSION=10.2.0
            % LMOD_FAMILY_MATLAB=matlab
            % LMOD_FAMILY_MATLAB_VERSION=R2024b
            % LMOD_FAMILY_SPM=spm
            % LMOD_FAMILY_SPM_VERSION=12
            % LMOD_FAMILY_WORKBENCH=workbench
            % LMOD_FAMILY_WORKBENCH_VERSION=2.0.1
            % LMOD_PKG=/usr/share/lmod/lmod
            % LMOD_ROOT=/usr/share/lmod
            % LMOD_SETTARG_FULL_SUPPORT=no
            % LMOD_VERSION=8.7.55
            % LMOD_sys=Linux
            % LOADEDMODULES=afni/25.1.03:ants/2.6.0:freesurfer/7.4.1:fsl/6.0.7.18:gcc/10.2.0:matlab/R2024b:spm/12:workbench/2.0.1
            % LOCAL=/home/jjlee/.local
            % LS_COLORS=
            % MANPATH=/usr/share/lmod/lmod/share/man:
            % MATLAB_MEM_MGR_VALUE=
            % MDCE_FORCE_MPI_OPTION=mpich
            % MDCE_PARALLEL=1
            % MEMKIND_HEAP_MANAGER=TBB
            % MKL_DOMAIN_NUM_THREADS=
            % MKL_NUM_THREADS=
            % MLM_WEB_LICENSE=false
            % MODULEPATH=/etc/modulefiles:/usr/share/modulefiles:/export/modulefiles:/opt/modulefiles
            % MODULEPATH_ROOT=/usr/share/modulefiles
            % MODULESHOME=/usr/share/lmod/lmod
            % MPIR_CVAR_CH3_INTERFACE_HOSTNAME=node16
            % MPI_LOCALNRANKS=2
            % MPI_LOCALRANKID=0
            % NILSRC=/home/jjlee/.local/lin64-nilsrc
            % OLDPWD=/matlab-local/R2024b
            % OMPI_MCA_orte_precondition_transports=0058c21f00000000-0058c21f00000000
            % OSG_LD_LIBRARY_PATH=/matlab-local/R2024b/sys/openscenegraph/lib/glnxa64
            % PARALLEL_SERVER_CMR=/matlab-local/R2024b
            % PARALLEL_SERVER_DEBUG=false
            % PARALLEL_SERVER_DECODE_FUNCTION=parallel.cluster.generic.communicatingDecodeFcn
            % PARALLEL_SERVER_JOB_LOCATION=Job13
            % PARALLEL_SERVER_MATLAB_ARGS=-parallel
            % PARALLEL_SERVER_MATLAB_EXE=/matlab-local/R2024b/bin/worker
            % PARALLEL_SERVER_NUM_THREADS=1
            % PARALLEL_SERVER_STORAGE_CONSTRUCTOR=makeFileStorageObject
            % PARALLEL_SERVER_STORAGE_LOCATION=/home/jjlee/.matlab/3p_cluster_jobs/chpc/twistor.attlocal.net.dhcp.wustl.edu/R2024b/nonshared
            % PARALLEL_SERVER_TOTAL_TASKS=2
            % PATH=/export/workbench:/usr/bin:/export/spm:/matlab-local/R2024b/bin:/export/gcc:/export/fsl:/export/freesurfer:/export/ants-2:/export/afni:/home/jjlee/julia-1.10.0/bin:/home/jjlee/miniconda/bin:/home/jjlee/.local/bin:/home/jjlee/.local/lin64-tools:/usr/local/bin
            % PMIX_BFROP_BUFFER_TYPE=PMIX_BFROP_BUFFER_NON_DESC
            % PMIX_DSTORE_21_BASE_PATH=/var/spool/slurmd/pmix.5816863.0//pmix_dstor_ds21_1199886
            % PMIX_DSTORE_ESH_BASE_PATH=/var/spool/slurmd/pmix.5816863.0//pmix_dstor_ds12_1199886
            % PMIX_GDS_MODULE=ds21,ds12,hash
            % PMIX_HOSTNAME=node16
            % PMIX_NAMESPACE=slurm.pmix.5816863.0
            % PMIX_RANK=0
            % PMIX_SECURITY_MODE=munge,native
            % PMIX_SERVER_TMPDIR=/var/spool/slurmd/pmix.5816863.0/
            % PMIX_SERVER_URI21=pmix-server.1199886;tcp4://127.0.0.1:58015
            % PMIX_SERVER_URI2=pmix-server.1199886;tcp4://127.0.0.1:58015
            % PMIX_SERVER_URI3=pmix-server.1199886;tcp4://127.0.0.1:58015
            % PMIX_SERVER_URI41=pmix-server.1199886;tcp4://127.0.0.1:58015
            % PMIX_SERVER_URI4=pmix-server.1199886;tcp4://127.0.0.1:58015
            % PMIX_SYSTEM_TMPDIR=/tmp
            % PMIX_VERSION=4.2.9
            % PMI_FD=6
            % PMI_RANK=0
            % PMI_SIZE=2
            % PRE_LD_PRELOAD=
            % PWD=/scratch/jjlee/Singularity/CCIR_01211
            % REFDIR=/home/jjlee/.local/atlas
            % RELEASE=/home/jjlee/.local/lin64-tools
            % SHARED=/scratch/jjlee/shared
            % SHELL=/bin/bash
            % SHLVL=1
            % SINGULARITY_HOME=/scratch/jjlee/Singularity
            % SLURMD_DEBUG=2
            % SLURMD_NODENAME=node16
            % SLURM_CLUSTER_NAME=chpc3
            % SLURM_CONF=/etc/slurm/slurm.conf
            % SLURM_CPUS_ON_NODE=2
            % SLURM_CPUS_PER_TASK=1
            % SLURM_CPU_BIND=quiet,mask_cpu:0x0C000000
            % SLURM_CPU_BIND_LIST=0x0C000000
            % SLURM_CPU_BIND_TYPE=mask_cpu:
            % SLURM_CPU_BIND_VERBOSE=quiet
            % SLURM_DISTRIBUTION=cyclic
            % SLURM_EXPORT_ENV=NONE
            % SLURM_GET_USER_ENV=1
            % SLURM_GTIDS=0
            % SLURM_JOBID=5816863
            % SLURM_JOB_ACCOUNT=manu_goyal
            % SLURM_JOB_CPUS_PER_NODE=2
            % SLURM_JOB_END_TIME=1753997383
            % SLURM_JOB_GID=1000070
            % SLURM_JOB_ID=5816863
            % SLURM_JOB_NAME=MATLAB_R2024b_Job13
            % SLURM_JOB_NODELIST=node16
            % SLURM_JOB_NUM_NODES=1
            % SLURM_JOB_PARTITION=tier1_cpu
            % SLURM_JOB_QOS=normal
            % SLURM_JOB_START_TIME=1753996783
            % SLURM_JOB_UID=1002720
            % SLURM_JOB_USER=jjlee
            % SLURM_LAUNCH_NODE_IPADDR=10.1.1.16
            % SLURM_LOCALID=0
            % SLURM_MEM_BIND=quiet,none
            % SLURM_MEM_BIND_LIST=
            % SLURM_MEM_BIND_TYPE=none
            % SLURM_MEM_BIND_VERBOSE=quiet
            % SLURM_MEM_PER_CPU=16384
            % SLURM_MPI_TYPE=pmix_v4
            % SLURM_NNODES=1
            % SLURM_NODEID=0
            % SLURM_NODELIST=node16
            % SLURM_NPROCS=1
            % SLURM_NTASKS=1
            % SLURM_NTASKS_PER_CORE=1
            % SLURM_OOM_KILL_STEP=0
            % SLURM_PMIXP_ABORT_AGENT_PORT=41381
            % SLURM_PMIX_MAPPING_SERV=(vector,(0,1,1))
            % SLURM_PRIO_PROCESS=0
            % SLURM_PROCID=0
            % SLURM_SRUN_COMM_HOST=10.1.1.16
            % SLURM_SRUN_COMM_PORT=46741
            % SLURM_STEPID=0
            % SLURM_STEP_ID=0
            % SLURM_STEP_LAUNCHER_PORT=46741
            % SLURM_STEP_NODELIST=node16
            % SLURM_STEP_NUM_NODES=1
            % SLURM_STEP_NUM_TASKS=1
            % SLURM_STEP_TASKS_PER_NODE=1
            % SLURM_SUBMIT_DIR=/home/jjlee
            % SLURM_SUBMIT_HOST=login01.cluster
            % SLURM_TASKS_PER_NODE=1
            % SLURM_TASK_PID=1199959
            % SLURM_TOPOLOGY_ADDR=node16
            % SLURM_TOPOLOGY_ADDR_PATTERN=node
            % SLURM_TRES_PER_TASK=cpu=1
            % SLURM_UMASK=0022
            % S_COLORS=auto
            % TMP=/scratch/jjlee/tmp
            % TMPDIR=/scratch/jjlee/tmp
            % TOOLBOX=/matlab-local/R2024b/toolbox
            % TZ=America/Chicago
            % USER=jjlee
            % XDG_DATA_DIRS=/home/jjlee/.local/share/flatpak/exports/share:/var/lib/flatpak/exports/share:/usr/local/share:/usr/share
            % XFILESEARCHPATH=/matlab-local/R2024b/sys/java/jre/glnxa64/jre/lib/locale/%L/%T/%N%S:
            % ZES_ENABLE_SYSMAN=1
            % _LMFILES_=/export/modulefiles/afni/25.1.03.lua:/export/modulefiles/ants/2.6.0.lua:/export/modulefiles/freesurfer/7.4.1.lua:/export/modulefiles/fsl/6.0.7.18.lua:/export/modulefiles/gcc/10.2.0.lua:/export/modulefiles/matlab/R2024b.lua:/export/modulefiles/spm/12.lua:/export/modulefiles/workbench/2.0.1.lua
            % _ModuleTable001_=X01vZHVsZVRhYmxlXyA9IHsKTVR2ZXJzaW9uID0gMywKY19yZWJ1aWxkVGltZSA9IGZhbHNlLApjX3Nob3J0VGltZSA9IGZhbHNlLApkZXB0aFQgPSB7fSwKZmFtaWx5ID0gewphZm5pID0gImFmbmkiLAphbnRzID0gImFudHMiLApmcmVlc3VyZmVyID0gImZyZWVzdXJmZXIiLApmc2wgPSAiZnNsIiwKZ2NjID0gImdjYyIsCm1hdGxhYiA9ICJtYXRsYWIiLApzcG0gPSAic3BtIiwKd29ya2JlbmNoID0gIndvcmtiZW5jaCIsCn0sCm1UID0gewphZm5pID0gewpmbiA9ICIvZXhwb3J0L21vZHVsZWZpbGVzL2FmbmkvMjUuMS4wMy5sdWEiLApmdWxsTmFtZSA9ICJhZm5pLzI1LjEuMDMiLApsb2FkT3JkZXIgPSAxLApwcm9wVCA9IHt9LApzdGFja0RlcHRoID0gMCwKc3RhdHVzID0g
            % _ModuleTable002_=ImFjdGl2ZSIsCnVzZXJOYW1lID0gImFmbmkiLAp3ViA9ICIwMDAwMDAwMjUuMDAwMDAwMDAxLjAwMDAwMDAwMy4qemZpbmFsIiwKfSwKYW50cyA9IHsKZm4gPSAiL2V4cG9ydC9tb2R1bGVmaWxlcy9hbnRzLzIuNi4wLmx1YSIsCmZ1bGxOYW1lID0gImFudHMvMi42LjAiLApsb2FkT3JkZXIgPSAyLApwcm9wVCA9IHt9LApzdGFja0RlcHRoID0gMCwKc3RhdHVzID0gImFjdGl2ZSIsCnVzZXJOYW1lID0gImFudHMiLAp3ViA9ICIwMDAwMDAwMDIuMDAwMDAwMDA2Lip6ZmluYWwiLAp9LApmcmVlc3VyZmVyID0gewpmbiA9ICIvZXhwb3J0L21vZHVsZWZpbGVzL2ZyZWVzdXJmZXIvNy40LjEubHVhIiwKZnVsbE5hbWUgPSAiZnJlZXN1cmZlci83LjQuMSIsCmxvYWRPcmRlciA9IDMs
            % _ModuleTable003_=CnByb3BUID0ge30sCnN0YWNrRGVwdGggPSAwLApzdGF0dXMgPSAiYWN0aXZlIiwKdXNlck5hbWUgPSAiZnJlZXN1cmZlciIsCndWID0gIjAwMDAwMDAwNy4wMDAwMDAwMDQuMDAwMDAwMDAxLip6ZmluYWwiLAp9LApmc2wgPSB7CmZuID0gIi9leHBvcnQvbW9kdWxlZmlsZXMvZnNsLzYuMC43LjE4Lmx1YSIsCmZ1bGxOYW1lID0gImZzbC82LjAuNy4xOCIsCmxvYWRPcmRlciA9IDQsCnByb3BUID0ge30sCnN0YWNrRGVwdGggPSAwLApzdGF0dXMgPSAiYWN0aXZlIiwKdXNlck5hbWUgPSAiZnNsIiwKd1YgPSAiMDAwMDAwMDA2LjAwMDAwMDAwMC4wMDAwMDAwMDcuMDAwMDAwMDE4Lip6ZmluYWwiLAp9LApnY2MgPSB7CmZuID0gIi9leHBvcnQvbW9kdWxlZmlsZXMvZ2NjLzEwLjIu
            % _ModuleTable004_=MC5sdWEiLApmdWxsTmFtZSA9ICJnY2MvMTAuMi4wIiwKbG9hZE9yZGVyID0gNSwKcHJvcFQgPSB7fSwKc3RhY2tEZXB0aCA9IDAsCnN0YXR1cyA9ICJhY3RpdmUiLAp1c2VyTmFtZSA9ICJnY2MiLAp3ViA9ICIwMDAwMDAwMTAuMDAwMDAwMDAyLip6ZmluYWwiLAp9LAptYXRsYWIgPSB7CmZuID0gIi9leHBvcnQvbW9kdWxlZmlsZXMvbWF0bGFiL1IyMDI0Yi5sdWEiLApmdWxsTmFtZSA9ICJtYXRsYWIvUjIwMjRiIiwKbG9hZE9yZGVyID0gNiwKcHJvcFQgPSB7fSwKc3RhY2tEZXB0aCA9IDAsCnN0YXR1cyA9ICJhY3RpdmUiLAp1c2VyTmFtZSA9ICJtYXRsYWIiLAp3ViA9ICIqci4wMDAwMDIwMjQuKmIuKnpmaW5hbCIsCn0sCnNwbSA9IHsKZm4gPSAiL2V4cG9ydC9tb2R1bGVm
            % _ModuleTable005_=aWxlcy9zcG0vMTIubHVhIiwKZnVsbE5hbWUgPSAic3BtLzEyIiwKbG9hZE9yZGVyID0gNywKcHJvcFQgPSB7fSwKc3RhY2tEZXB0aCA9IDAsCnN0YXR1cyA9ICJhY3RpdmUiLAp1c2VyTmFtZSA9ICJzcG0iLAp3ViA9ICIwMDAwMDAwMTIuKnpmaW5hbCIsCn0sCndvcmtiZW5jaCA9IHsKZm4gPSAiL2V4cG9ydC9tb2R1bGVmaWxlcy93b3JrYmVuY2gvMi4wLjEubHVhIiwKZnVsbE5hbWUgPSAid29ya2JlbmNoLzIuMC4xIiwKbG9hZE9yZGVyID0gOCwKcHJvcFQgPSB7fSwKc3RhY2tEZXB0aCA9IDAsCnN0YXR1cyA9ICJhY3RpdmUiLAp1c2VyTmFtZSA9ICJ3b3JrYmVuY2giLAp3ViA9ICIwMDAwMDAwMDIuMDAwMDAwMDAwLjAwMDAwMDAwMS4qemZpbmFsIiwKfSwKfSwKbXBhdGhB
            % _ModuleTable006_=ID0gewoiL2V0Yy9tb2R1bGVmaWxlcyIKLCAiL3Vzci9zaGFyZS9tb2R1bGVmaWxlcyIKLCAiL2V4cG9ydC9tb2R1bGVmaWxlcyIsICIvb3B0L21vZHVsZWZpbGVzIiwKfSwKc3lzdGVtQmFzZU1QQVRIID0gIi9ldGMvbW9kdWxlZmlsZXM6L3Vzci9zaGFyZS9tb2R1bGVmaWxlczovZXhwb3J0L21vZHVsZWZpbGVzOi9vcHQvbW9kdWxlZmlsZXMiLAp9Cg==
            % _ModuleTable_Sz_=6
            % __LMOD_REF_COUNT_MODULEPATH=/etc/modulefiles:1;/usr/share/modulefiles:1;/export/modulefiles:1;/opt/modulefiles:1
            % __LMOD_REF_COUNT_PATH=/export/workbench:1;/usr/bin:9;/export/spm:1;/matlab-local/R2024b/bin:1;/export/gcc:1;/export/fsl:1;/export/freesurfer:1;/export/ants-2:1;/export/afni:1;/home/jjlee/julia-1.10.0/bin:1;/home/jjlee/miniconda/bin:1;/home/jjlee/.local/bin:1;/home/jjlee/.local/lin64-tools:1;/usr/local/bin:1
            % which_declare=declare -f
            % }
            % }
            % }

            arguments
                globbing_mat {mustBeFile} = ...
                    fullfile( ...
                    getenv("HOME"), ...
                    "mnt", "CHPC_scratch", "Singularity", "CCIR_01211", "mlvg_Lee2025Par_globbing_fdg.mat")
                opts.globbing_var = "globbed"
                opts.selection_indices double = 2  % total ~ 1:58 for ho, 1:69 for co, 1:112 for oo, 1:57 for fdg
                opts.account {mustBeTextScalar} = "manu_goyal"
                opts.M {mustBeInteger} = 1
            end

            ld = load(globbing_mat);
            globbed = convertCharsToStrings(ld.(opts.globbing_var));
            if ~isempty(opts.selection_indices)
                globbed = globbed(opts.selection_indices);
            end
            globbed = strrep(globbed, "createNiftiMovingAvgFrames", "createNiftiStatic");
            fprintf("%s:globbed:\n", stackstr())
            disp(size(globbed))
            disp(ascol(globbed))

            % contact cluster slurm

            c = mlvg.CHPC3.propcluster(opts.account, mempercpu='16gb', walltime='0:10:00');
            disp(c)
            disp(c.AdditionalProperties)
            try
                j = c.batch( ...
                    @mlvg.Lee2025Par.test_fsl, ...
                    1, ...
                    {globbed}, ...
                    'Pool', opts.M, ...
                    'CurrentFolder', '/scratch/jjlee/Singularity/CCIR_01211', ...
                    'AutoAddClientPath', false);
            catch ME
                handwarning(ME)
            end

            [msg,id] = lastwarn();
        end

        function [j,c,msg,id] = cluster_build_schaeffer_parc(globbing_mat, opts)
            %% for clusters running Matlab parallel server

            arguments
                globbing_mat {mustBeFile} = ...
                    fullfile( ...
                    getenv("SINGULARITY_HOME"), "CCIR_01211", "mlvg_Lee2025Par_globbing_co.mat")
                opts.globbing_var = "globbed"
                opts.selection_indices double = []  % total ~ 1:58 for ho, 1:69 for co, 1:112 for oo
                opts.Ncol {mustBeInteger} = 4
                opts.account {mustBeTextScalar} = "manu_goyal"
            end
            ld = load(globbing_mat);
            globbed = convertCharsToStrings(ld.(opts.globbing_var));
            globbed = asrow(globbed);
            if ~isempty(opts.selection_indices)
                globbed = globbed(opts.selection_indices);
            end
            is_delayed = contains(globbed(1), "-delay");
            if contains(globbing_mat, "co", IgnoreCase=true)
                c = mlvg.CHPC3.propcluster(opts.account, mempercpu='128gb', walltime='12:00:00');
                opts.Ncol = 2;
            elseif contains(globbing_mat, "fdg", IgnoreCase=true)
                if ~is_delayed
                    c = mlvg.CHPC3.propcluster(opts.account, mempercpu='128gb', walltime='12:00:00');
                    opts.Ncol = 2;
                else
                    c = mlvg.CHPC3.propcluster(opts.account, mempercpu='16gb', walltime='2:00:00');
                end
            else
                if ~is_delayed
                    c = mlvg.CHPC3.propcluster(opts.account, mempercpu='64gb', walltime='12:00:00');
                else
                    c = mlvg.CHPC3.propcluster(opts.account, mempercpu='16gb', walltime='1:00:00');
                end
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

            disp(c.AdditionalProperties)
            for irow = 1:Nrow
                try
                    j = c.batch( ...
                        @mlvg.Lee2025Par.par_build_schaeffer_parc, ...
                        1, ...
                        {globbed(irow, :)}, ...
                        'Pool', opts.Ncol, ...
                        'CurrentFolder', '/scratch/jjlee/Singularity/CCIR_01211', ...
                        'AutoAddClientPath', false);
                catch ME
                    handwarning(ME)
                end
            end

            [msg,id] = lastwarn();
        end

        function [j,c,msg,id] = cluster_call(globbing_mat, opts)
            %% for clusters running Matlab parallel server

            arguments
                globbing_mat {mustBeFile} = ...
                    fullfile( ...
                    getenv("SINGULARITY_HOME"), "CCIR_01211", "srcdata_fdg.mat")
                opts.globbing_var = "srcdata_fdg"
                opts.selection_indices double = [1, 3:78]  % total ~ 1:58 for ho, 1:69 for co, 1:112 for oo
                opts.Ncol {mustBeInteger} = 4
                opts.method {mustBeTextScalar} = "do_make_input_func"
                opts.reference_tracer {mustBeTextScalar} = "fdg"
                opts.steps {mustBeNumericOrLogical} = 1
                opts.account {mustBeTextScalar} = "manu_goyal"
            end
            c = mlvg.CHPC3.propcluster(opts.account, mempercpu='128gb', walltime='3:00:00');
            disp(c.AdditionalProperties)
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

        function [j,c,msg,id] = cluster_time_align_static(globbing_mat, opts)
            %% for clusters running Matlab parallel server

            arguments
                globbing_mat {mustBeFile} = ...
                    fullfile( ...
                    getenv("HOME"), ...
                    "mnt", "CHPC_scratch", "Singularity", "CCIR_01211", "mlvg_Lee2025Par_globbing_fdg.mat")
                opts.globbing_var = "globbed"
                opts.selection_indices double = 1:4  % total ~ 1:58 for ho, 1:69 for co, 1:112 for oo, 1:57 for fdg
                opts.Ncol {mustBeInteger} = 4
                opts.account {mustBeTextScalar} = "manu_goyal"
            end
            ld = load(globbing_mat);
            globbed = convertCharsToStrings(ld.(opts.globbing_var));
            globbed = asrow(globbed);
            if ~isempty(opts.selection_indices)
                globbed = globbed(opts.selection_indices);
            end
            if contains(globbing_mat, "fdg", IgnoreCase=true) || ...
                    contains(globbing_mat, "co", IgnoreCase=true)
                c = mlvg.CHPC3.propcluster(opts.account, mempercpu='48gb', walltime='4:00:00');
                opts.Ncol = 4;
            else
                c = mlvg.CHPC3.propcluster(opts.account, mempercpu='32gb', walltime='4:00:00');
                opts.Ncol = 4;
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

            disp(c.AdditionalProperties)
            for irow = 1:Nrow
                try
                    j = c.batch( ...
                        @mlvg.Lee2025Par.par_time_align_static, ...
                        1, ...
                        {globbed(irow, :), "M", opts.Ncol}, ...
                        'Pool', opts.Ncol, ...
                        'CurrentFolder', '/scratch/jjlee/Singularity/CCIR_01211', ...
                        'AutoAddClientPath', false);
                catch ME
                    handwarning(ME)
                end
            end

            [msg,id] = lastwarn();
        end

        function [j,c,msg,id] = cluster_reflirt_t1w(globbing_mat, opts)
            %% for clusters running Matlab parallel server

            arguments
                globbing_mat {mustBeFile} = ...
                    fullfile(getenv("HOME"), "mnt", "CHPC_scratch", "Singularity", "CCIR_01211", "srcdata_ho_todo.mat")
                opts.globbing_var = "srcdata_todo"
                opts.selection_indices double = []
                opts.Ncol {mustBeInteger} = 4
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
            c = mlvg.CHPC3.propcluster(opts.account, mempercpu='128gb', walltime='4:00:00');
            disp(c.AdditionalProperties)
            for irow = 1:Nrow
                try
                    j = c.batch( ...
                        @mlvg.Lee2025Par.par_reflirt_t1w, ...
                        1, ...
                        {globbed(irow, :)}, ...
                        'Pool', opts.Ncol, ...
                        'CurrentFolder', '/scratch/jjlee/Singularity/CCIR_01211', ...
                        'AutoAddClientPath', false);
                catch ME
                    handwarning(ME)
                end
            end

            [msg,id] = lastwarn();
        end

        function durations = par_build_schaeffer_parc(fqfns, opts)
            %% e.g.,
            %  sub-108306_ses-20230227113853_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz ->
            %  sub-108306_ses-20230227113853_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames-ParcSchaeffer-reshape-to-schaeffer-schaeffer.nii.gz
            %  for OO, uses 70 GB of memory per parallel process

            arguments
                fqfns {mustBeText}
                opts.out_dir {mustBeFolder} = "/scratch/jjlee/Singularity/CCIR_01211"
                opts.M {mustBeScalarOrEmpty} = 4  % 8
                opts.select {mustBeInteger} = 1
                opts.do_make_finite logical = false
            end
            if ~contains(fqfns(1), opts.out_dir)
                fqfns = fullfile(opts.out_dir, fqfns);
            end

            import mlkinetics.*

            durations = nan(1, length(fqfns));
            
            if isscalar(fqfns)
                try
                    % setup
                    mlvg.CHPC3.setenvs();

                    fqfn = fqfns(1);
                    durations = mlvg.Lee2025Par.build_schaeffer_parc( ...
                        fqfn, ...
                        out_dir=opts.out_dir, do_make_finite=opts.do_make_finite);
                catch ME
                    handwarning(ME)
                end
            else
                parfor (fidx = 1:length(fqfns), opts.M)
                    try
                        % setup
                        mlvg.CHPC3.setenvs();

                        fqfn = fqfns(fidx);
                        durations(fidx) = mlvg.Lee2025Par.build_schaeffer_parc( ...
                            fqfn, ...
                            out_dir=opts.out_dir, do_make_finite=opts.do_make_finite); %#ok<PFBNS>
                    catch ME
                        handwarning(ME)
                    end
                end
            end
        end

        function durations = par_call_ifk(nii, opts)
            arguments
                nii {mustBeText} = "sourcedata/sub-108007/ses-20210219145054/pet/sub-108007_ses-20210219145054_trc-ho_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz"
                opts.out_dir {mustBeFolder} = "/scratch/jjlee/Singularity/CCIR_01211"
                opts.method {mustBeTextScalar} = "do_make_input_func"
                opts.reference_tracer {mustBeTextScalar} = "fdg"
                opts.steps {mustBeNumericOrLogical} = 1
            end
            
            nii = nii(~arrayfun(@isempty, nii));  % Remove empty cells
            durations = nan(1, length(nii));

            parfor sidx = 1:length(nii)

                tic;
            
                % setup
                mlvg.CHPC3.setenvs();

                try
                    % exclusions
                    if 5 == opts.steps %#ok<PFBNS>
                        fn = extractBefore(mybasename(nii(sidx)), "_proc") + "_proc-MipIdif_idif.nii.gz";
                        target = fullfile(fileparts(nii(sidx)), fn);
                        target = strrep(target, "sourcedata", "derivatives");
                        if isfile(target)
                            fprintf("%s: skipping existing %s\n", stackstr(), target);
                            durations(sidx) = toc;
                            continue
                        end
                    end

                    % construct & call
                    lp = mlvg.Lee2025Par(nii(sidx), out_dir=opts.out_dir); %#ok<PFBNS>
                    call_ifk(lp, method=opts.method, steps=opts.steps, reference_tracer=opts.reference_tracer);
                catch ME
                    handwarning(ME)
                end

                durations(sidx) = toc;
            end
        end

        function durations = par_time_align_static(nii, opts)
            arguments
                nii {mustBeText} = "sourcedata/sub-108007/ses-20210219143132/pet/sub-108007_ses-20210219143132_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz"
                opts.out_dir {mustBeFolder} = "/scratch/jjlee/Singularity/CCIR_01211"
                opts.M {mustBeInteger} = 8
            end
            
            nii = nii(~arrayfun(@isempty, nii));  % Remove empty cells
            durations = nan(1, length(nii));
            
            if 1 == opts.M
                for sidx = 1:length(nii)

                    tic;

                    % setup
                    mlvg.CHPC3.setenvs();

                    try
                        nii_fqfn = fullfile(opts.out_dir, nii(sidx)); 
                        globbed = mglob(strrep(nii_fqfn, "delay0", "delay*"));
                        mlvg.Lee2025.time_align_static(globbed);
                    catch ME
                        handwarning(ME)
                    end

                    durations(sidx) = toc;
                end
                return
            end

            parfor (sidx = 1:length(nii), opts.M)

                tic;
            
                % setup
                mlvg.CHPC3.setenvs();

                try
                    nii_fqfn = fullfile(opts.out_dir, nii(sidx)); %#ok<PFBNS>
                    globbed = mglob(strrep(nii_fqfn, "delay0", "delay*"));
                    mlvg.Lee2025.time_align(globbed);
                catch ME
                    handwarning(ME)
                end

                durations(sidx) = toc;
            end
        end

        function durations = par_reflirt_t1w(nii, opts)
            arguments
                nii {mustBeText} = "sourcedata/sub-108007/ses-20210219143132/pet/sub-108007_ses-20210219143132_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz"
                opts.out_dir {mustBeFolder} = "/scratch/jjlee/Singularity/CCIR_01211"
                opts.M {mustBeNumeric} = []
            end
            if isempty(opts.M)
                opts.M = length(nii);
            end
            
            nii = nii(~arrayfun(@isempty, nii));  % Remove empty cells
            durations = nan(1, length(nii));

            parfor (sidx = 1:length(nii), opts.M)

                tic;
            
                % setup
                mlvg.CHPC3.setenvs();

                try
                    nii_fqfn = fullfile(opts.out_dir, nii(sidx)); %#ok<PFBNS>
                    mlvg.Lee2025.flirt_t1w(nii_fqfn); 
                catch ME
                    handwarning(ME)
                end

                durations(sidx) = toc;
            end
        end

        function durations = par_align_delayed_static(nii, opts)
            arguments
                nii {mustBeText} = "sourcedata/sub-108007/ses-20210219143132/pet/sub-108007_ses-20210219143132_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz"
                opts.out_dir {mustBeFolder} = "/scratch/jjlee/Singularity/CCIR_01211"
            end
            
            nii = nii(~arrayfun(@isempty, nii));  % Remove empty cells
            durations = nan(1, length(nii));

            parfor sidx = 1:length(nii)

                tic;
            
                % setup
                mlvg.CHPC3.setenvs();

                try
                    nii_fqfn = fullfile(opts.out_dir, nii(sidx)); %#ok<PFBNS>
                    mlvg.Lee2025.flirt_t1w(nii_fqfn); 
                catch ME
                    handwarning(ME)
                end

                durations(sidx) = toc;
            end
        end
    
        function r = test_fsl(nii, opts)
            arguments
                nii {mustBeText} = "sourcedata/sub-108007/ses-20210219143132/pet/sub-108007_ses-20210219143132_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames.nii.gz"
                opts.out_dir {mustBeFolder} = "/scratch/jjlee/Singularity/CCIR_01211"
                opts.M {mustBeInteger} = 1
            end

            r = "";

            if 1 == opts.M
                for sidx = 1:length(nii)

                    % setup
                    mlvg.CHPC3.setenvs();

                    try
                        %[~, result] = system('env | sort');
                        %[~, result] = system('~/bin/mlenv.sh env | sort');
                        %[~, result] = mysystem('env | sort');
                        %[~,result] = mysystem("flirt");
                        [~,result] = mysystem("fslhd " + nii(sidx));
                    catch ME
                        handwarning(ME)
                    end

                    r(sidx) = result;
                end
                return
            end

            parfor (sidx = 1:length(nii), opts.M)

                % setup
                mlvg.CHPC3.setenvs();

                try
                    result = mysystem("fslhd " + nii(sidx));
                catch ME
                    handwarning(ME)
                end

                r(sidx) = result;
            end

        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
