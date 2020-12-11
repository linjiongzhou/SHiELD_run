#!/bin/tcsh -f
#SBATCH --output=/lustre/f2/scratch/Linjiong.Zhou/SHiELD/stdout/%x.o%j
#SBATCH --job-name=C48_20150801.00Z
#SBATCH --partition=batch
#SBATCH --account=gfdl_w
#SBATCH --time=02:00:00
#SBATCH --cluster=c3
#SBATCH --nodes=6
#SBATCH --export=NAME=20150801.00Z,MEMO=_RT2021,EXE=x,ALL

# This script is optimized for GFDL MP runs using GFS ICs
# Linjiong.Zhou@noaa.gov

set echo

set BASEDIR    = "/lustre/f2/scratch/${USER}/SHiELD"
set INPUT_DATA = "/lustre/f2/pdata/gfdl/gfdl_W/fvGFS_INPUT_DATA"
set BUILD_AREA = "/lustre/f2/dev/${USER}/SHiELD/SHiELD_build"
set RUN_AREA = "/lustre/f2/dev/${USER}/SHiELD/SHiELD_run"

# release number for the script
set RELEASE = "`cat ${BUILD_AREA}/release`"

# case specific details
set TYPE = "nh"         # choices:  nh, hydro
set MODE = "32bit"      # choices:  32bit, 64bit
set MONO = "non-mono"   # choices:  mono, non-mono
set CASE = "C48"
#set NAME = "20150801.00Z"
#set MEMO = "_RT2018"
#set EXE = "x"
set HYPT = "on"         # choices:  on, off  (controls hyperthreading)
set COMP = "prod"       # choices:  debug, repro, prod
set NO_SEND = "send"    # choices:  send, no_send

# directory structure
set WORKDIR    = ${BASEDIR}/${RELEASE}/${NAME}.${CASE}.${TYPE}.${MODE}.${MONO}${MEMO}/
set executable = ${BUILD_AREA}/Build/bin/SHiELD_${TYPE}.${COMP}.${MODE}.${EXE}

# input filesets
set ICS  = ${INPUT_DATA}/global.v202012/${CASE}/${NAME}_IC
set FIX  = ${INPUT_DATA}/fix.v201912
set GRID = ${INPUT_DATA}/global.v202012/${CASE}/GRID
set FIX_bqx  = ${INPUT_DATA}/climo_data.v201807
set FIX_sfc = ${GRID}/fix_sfc

# sending file to gfdl
set gfdl_archive = /archive/${USER}/SHiELD/${RELEASE}/${NAME}.${CASE}.${TYPE}.${MODE}.${MONO}${MEMO}/
set SEND_FILE =  ${BUILD_AREA}/site/send_file_slurm.csh
set TIME_STAMP = ${BUILD_AREA}/site/time_stamp.csh

# changeable parameters
    # dycore definitions
    set npx = "49"
    set npy = "49"
    set npz = "91"
    set layout_x = "4"
    set layout_y = "4" 
    set io_layout = "1,1"
    set nthreads = "4"

    # blocking factor used for threading and general physics performance
    set blocksize = "32"

    # run length
    set months = "0"
    set days = "10"
    set hours = "0"
    set dt_atmos = "150"

    # set the pre-conditioning of the solution
    # =0 implies no pre-conditioning
    # >0 means new adiabatic pre-conditioning
    # <0 means older adiabatic pre-conditioning
    set na_init = 0

    # variables for controlling initialization of NCEP/NGGPS ICs
    set filtered_terrain = ".true."
    set ncep_levs = "64"
    set gfs_dwinds = ".true."

    # variables for gfs diagnostic output intervals and time to zero out time-accumulated data
#    set fdiag = "6.,12.,18.,24.,30.,36.,42.,48.,54.,60.,66.,72.,78.,84.,90.,96.,102.,108.,114.,120.,126.,132.,138.,144.,150.,156.,162.,168.,174.,180.,186.,192.,198.,204.,210.,216.,222.,228.,234.,240."
    set fdiag = "6."
    set fhzer = "6."
    set fhcyc = "24."

    # determines whether FV3 or GFS physics calculate geopotential
    set gfs_phil = ".false."

    # determine whether ozone production occurs in GFS physics
    set ozcalc = ".true."

    # set various debug options
    set no_dycore = ".false."
    set dycore_only = ".false."
    set chksum_debug = ".false."
    set print_freq = "6"

    if (${TYPE} == "nh") then
      # non-hydrostatic options
      set make_nh = ".F."
      set hydrostatic = ".F."
      set phys_hydrostatic = ".F."     # can be tested
      set use_hydro_pressure = ".F."   # can be tested
      set consv_te = "1."
        # time step parameters in FV3
      set k_split = "1"
      set n_split = "2"
    else
      # hydrostatic options
      set make_nh = ".F."
      set hydrostatic = ".T."
      set phys_hydrostatic = ".F."     # will be ignored in hydro mode
      set use_hydro_pressure = ".T."   # have to be .T. in hydro mode
      set consv_te = "0."
        # time step parameters in FV3
      set k_split = "2"
      set n_split = "6"
    endif

    if (${MONO} == "mono" || ${MONO} == "monotonic") then
      # monotonic options
      set d_con = "1."
      set do_vort_damp = ".false."
    else
      # non-monotonic options
      set d_con = "1."
      set do_vort_damp = ".true."
    endif

    # variables for hyperthreading
    if (${HYPT} == "on") then
      set hyperthread = ".true."
      set j_opt = "-j2"
	  set div = 2
    else
      set hyperthread = ".false."
      set j_opt = "-j1"
	  set div = 1
    endif

# when running with threads, need to use the following command
    @ npes = ${layout_x} * ${layout_y} * 6
	@ skip = ${nthreads} / ${div}
	set run_cmd = "srun --ntasks=$npes --cpus-per-task=$skip ./$executable:t"

	setenv SLURM_CPU_BIND verbose

    setenv MPICH_ENV_DISPLAY
    setenv MPICH_MPIIO_CB_ALIGN 2
    setenv MALLOC_MMAP_MAX_ 0
    setenv MALLOC_TRIM_THRESHOLD_ 536870912
    setenv NC_BLKSZ 1M
# necessary for OpenMP when using Intel
    setenv KMP_STACKSIZE 256m

\rm -rf $WORKDIR/rundir

mkdir -p $WORKDIR/rundir
cd $WORKDIR/rundir

mkdir -p RESTART

# build the date for curr_date and diag_table from NAME
unset echo
set y = `echo ${NAME} | cut -c1-4`
set m = `echo ${NAME} | cut -c5-6`
set d = `echo ${NAME} | cut -c7-8`
set h = `echo ${NAME} | cut -c10-11`
set echo
set curr_date = "${y},${m},${d},${h},0,0"

# build the diag_table with the experiment name and date stamp
cat > diag_table << EOF
${NAME}.${CASE}.${MODE}.${MONO}
$y $m $d $h 0 0 
EOF
cat ${RUN_AREA}/diag_table_6species >> diag_table

# copy over the other tables and executable
cp ${RUN_AREA}/data_table data_table
cp ${RUN_AREA}/field_table_6species field_table
cp $executable .

mkdir -p INPUT

# Grid and orography data
ln -sf ${GRID}/* INPUT/

# Date specific ICs
ln -sf ${ICS}/* INPUT/

# GFS FIX data
ln -sf $FIX/ozprdlos_2015_new_sbuvO3_tclm15_nuchem.f77 INPUT/global_o3prdlos.f77
ln -sf $FIX/global_h2o_pltc.f77 INPUT/global_h2oprdlos.f77
ln -sf $FIX/global_solarconstant_noaa_an.txt INPUT/solarconstant_noaa_an.txt
ln -sf $FIX/global_sfc_emissivity_idx.txt INPUT/sfc_emissivity_idx.txt
ln -sf $FIX/global_co2historicaldata_glob.txt INPUT/co2historicaldata_glob.txt
ln -sf $FIX/co2monthlycyc.txt INPUT/co2monthlycyc.txt
foreach file ( $FIX/fix_co2_proj/global_co2historicaldata_????.txt )
	ln -sf $file INPUT/`echo $file:t | sed s/global_co2historicaldata/co2historicaldata/g`
end
ln -sf $FIX/global_climaeropac_global.txt INPUT/aerosol.dat
foreach file ( $FIX/global_volcanic_aerosols_????-????.txt )
	ln -sf $file INPUT/`echo $file:t | sed s/global_volcanic_aerosols/volcanic_aerosols/g`
end

cat > input.nml <<EOF
 &amip_interp_nml
     interp_oi_sst = .true.
     use_ncep_sst = .true.
     use_ncep_ice = .false.
     no_anom_sst = .false.
     data_set = 'reynolds_oi',
     date_out_of_range = 'climo',
/

 &atmos_model_nml
     blocksize = $blocksize
     chksum_debug = $chksum_debug
     dycore_only = $dycore_only
     fdiag = $fdiag
     first_time_step = .false.
/

 &fms_io_nml
       checksum_required   = .false.
       max_files_r = 100,
       max_files_w = 100,
/

 &fms_nml
       clock_grain = 'ROUTINE',
       domains_stack_size = 3000000,
       print_memory_usage = .false.
/

 &fv_grid_nml
       grid_file = 'INPUT/grid_spec.nc'
/

 &fv_core_nml
       layout   = $layout_x,$layout_y
       io_layout = $io_layout
       npx      = $npx
       npy      = $npy
       ntiles   = 6
       npz    = $npz
       grid_type = -1
       make_nh = $make_nh
       fv_debug = .F.
       range_warn = .T.
       reset_eta = .F.
       n_sponge = 30
       nudge_qv = .T.
       rf_fast = .F.
       tau = 5.
       rf_cutoff = 7.5e2
       d2_bg_k1 = 0.15
       d2_bg_k2 = 0.02
       kord_tm = -9
       kord_mt =  9
       kord_wz =  9
       kord_tr =  9
       hydrostatic = $hydrostatic
       phys_hydrostatic = $phys_hydrostatic
       use_hydro_pressure = $use_hydro_pressure
       beta = 0.
       a_imp = 1.
       p_fac = 0.1
       k_split  = $k_split
       n_split  = $n_split
       nwat = 6 
       na_init = $na_init
       d_ext = 0.0
       dnats = 1
       fv_sg_adj = 600
       d2_bg = 0.
       nord =  3
       dddmp = 0.2
       d4_bg = 0.15
       vtdm4 = 0.03
       delt_max = 0.002
       ke_bg = 0.
       do_vort_damp = $do_vort_damp
       external_ic = .T.
       gfs_phil = $gfs_phil
       nggps_ic = .T.
       mountain = .F.
       ncep_ic = .F.
       d_con = $d_con
       hord_mt = 5
       hord_vt = 5
       hord_tm = 5
       hord_dp = -5
       hord_tr = -5
       adjust_dry_mass = .F.
       consv_te = $consv_te
       do_sat_adj = .F.
       consv_am = .F.
       fill = .T.
       dwind_2d = .F.
       print_freq = $print_freq
       warm_start = .F.
       no_dycore = $no_dycore
       z_tracer = .T.
       do_inline_mp = .T.
/

 &coupler_nml
       months = $months
       days  = $days
       hours = $hours
       dt_atmos = $dt_atmos
       dt_ocean = $dt_atmos
       current_date =  $curr_date
       calendar = 'julian'
       memuse_verbose = .false.
       atmos_nthreads = $nthreads
       use_hyper_thread = $hyperthread
/

 &external_ic_nml 
       filtered_terrain = $filtered_terrain
       levp = $ncep_levs
       gfs_dwinds = $gfs_dwinds
       checker_tr = .F.
       nt_checker = 0
/

 &gfs_physics_nml
       fhzero         = $fhzer
       ldiag3d        = .false.
       fhcyc          = $fhcyc
       nst_anl        = .true.
       use_ufo        = .true.
       pre_rad        = .false.
       ncld           = 5
       zhao_mic       = .false.
       pdfcld         = .true.
       fhswr          = 3600.
       fhlwr          = 3600.
       ialb           = 1
       iems           = 1
       IAER           = 111
       ico2           = 2
       isubc_sw       = 2
       isubc_lw       = 2
       isol           = 2
       lwhtr          = .true.
       swhtr          = .true.
       cnvgwd         = .true.
       do_deep        = .true.
       shal_cnv       = .true.
       cal_pre        = .false.
       redrag         = .true.
       dspheat        = .true.
       hybedmf        = .false.
       random_clds    = .false.
       trans_trac     = .true.
       cnvcld         = .false.
       imfshalcnv     = 2
       imfdeepcnv     = 2
       cdmbgwd        = 3.5, 0.25
       prslrd0        = 0.
       ivegsrc        = 1
       isot           = 1
       ysupbl         = .false.
       satmedmf       = .true.
       isatmedmf      = 1
       do_dk_hb19     = .false.
       xkzminv        = 0.0
	   xkzm_m         = 1.5
       xkzm_h         = 1.5
	   xkzm_ml        = 1.0
       xkzm_hl        = 1.0
	   xkzm_mi        = 1.5
       xkzm_hi        = 1.5
       cap_k0_land    = .false.
       cloud_gfdl     = .true.
       do_inline_mp   = .true.
       do_ocean       = .true.
       do_z0_hwrf17_hwonly = .true.
/

 &ocean_nml
     mld_option       = "obs"
     ocean_option     = "MLM"
     restore_method   = 2
     mld_obs_ratio    = 1.
     use_rain_flux    = .true.
     sst_restore_tscale = 2.
     start_lat        = -30.
     end_lat          = 30.
     Gam              = 0.2
     use_old_mlm      = .true.
     do_mld_restore   = .true.
	 mld_restore_tscale = 2.
     stress_ratio     = 1.
     eps_day          = 10.
/

 &gfdl_mp_nml
       sedi_transport = .true.
       do_sedi_w = .true.
       do_sedi_heat = .false.
       disp_heat = .true.
       rad_snow = .true.
       rad_graupel = .true.
       rad_rain = .true.
       const_vi = .false.
       const_vs = .false.
       const_vg = .false.
       const_vr = .false.
       vi_fac = 1.
       vs_fac = 1.
       vg_fac = 1.
       vr_fac = 1.
       vi_max = 1.
       vs_max = 2.
       vg_max = 12.
       vr_max = 12.
       qi_lim = 1.
       prog_ccn = .false.
       do_qa = .true.
       do_sat_adj = .false.
       tau_l2v = 225.
       tau_v2l = 150.
       tau_g2v = 900.
       rthresh = 10.e-6
       dw_land = 0.16
       dw_ocean = 0.10
       ql_gen = 1.0e-3
       ql_mlt = 1.0e-3
       qi0_crt = 8.0e-5
       qs0_crt = 1.0e-3
       tau_i2s = 1000.
       c_psaci = 0.05
       c_pgacs = 0.01
       rh_inc = 0.30
       rh_inr = 0.30
       rh_ins = 0.30
       ccn_l = 300.
       ccn_o = 200.
       c_paut = 0.5
       c_cracw = 0.8
       use_ppm = .false.
       mono_prof = .true.
       z_slope_liq = .true.
       z_slope_ice = .true.
       fix_negative = .true.
       icloud_f = 0
       do_cld_adj = .true.
       f_dq_p = 3.0
/

 &cld_eff_rad_nml
       qmin = 1.0e-12
       beta = 1.22
       rewflag = 1
       reiflag = 5
       rewmin = 5.0
       rewmax = 10.0
       reimin = 10.0
       reimax = 150.0
       rermin = 10.0
       rermax = 10000.0
       resmin = 150.0
       resmax = 10000.0
       liq_ice_combine = .false.
/

 &diag_manager_nml 
       prepend_date = .F.
/

  &interpolator_nml
       interp_method = 'conserve_great_circle'
/

&namsfc
       FNGLAC   = "$FIX/global_glacier.2x2.grb",
       FNMXIC   = "$FIX/global_maxice.2x2.grb",
       FNTSFC   = "$FIX/RTGSST.1982.2012.monthly.clim.grb",
       FNMLDC   = "$FIX_bqx/mld/mld_DR003_c1m_reg2.0.grb"
       FNSNOC   = "$FIX/global_snoclim.1.875.grb",
       FNZORC   = "igbp",
       FNALBC   = "$FIX_sfc/${CASE}.snowfree_albedo.tileX.nc",
       FNALBC2  = "$FIX_sfc/${CASE}.facsf.tileX.nc",
       FNAISC   = "$FIX/CFSR.SEAICE.1982.2012.monthly.clim.grb",
       FNTG3C   = "$FIX_sfc/${CASE}.substrate_temperature.tileX.nc",
       FNVEGC   = "$FIX_sfc/${CASE}.vegetation_greenness.tileX.nc",
       FNVETC   = "$FIX_sfc/${CASE}.vegetation_type.tileX.nc",
       FNSOTC   = "$FIX_sfc/${CASE}.soil_type.tileX.nc",
       FNSMCC   = "$FIX/global_soilmgldas.t1534.3072.1536.grb",
       FNMSKH   = "$FIX/global_slmask.t1534.3072.1536.grb",
       FNTSFA   = "",
       FNACNA   = "",
       FNSNOA   = "",
       FNVMNC   = "$FIX_sfc/${CASE}.vegetation_greenness.tileX.nc",
       FNVMXC   = "$FIX_sfc/${CASE}.vegetation_greenness.tileX.nc",
       FNSLPC   = "$FIX_sfc/${CASE}.slope_type.tileX.nc",
       FNABSC   = "$FIX_sfc/${CASE}.maximum_snow_albedo.tileX.nc",
       LDEBUG   =.false.,
       FSMCL(2) = 99999
       FSMCL(3) = 99999
       FSMCL(4) = 99999
       FTSFS    = 90
       FAISS    = 99999
       FSNOL    = 99999
       FSICL    = 99999
       FTSFL    = 99999,
       FAISL    = 99999,
       FVETL    = 99999,
       FSOTL    = 99999,
       FvmnL    = 99999,
       FvmxL    = 99999,
       FSLPL    = 99999,
       FABSL    = 99999,
       FSNOS    = 99999,
       FSICS    = 99999,
/
EOF

# run the executable
   ${run_cmd} | tee fms.out
   if ( $? != 0 ) then
     exit
   endif

if ($NO_SEND == "no_send") then
  exit
endif

#########################################################################
# generate date for file names
########################################################################

    set begindate = `$TIME_STAMP -bhf digital`
    if ( $begindate == "" ) set begindate = tmp`date '+%j%H%M%S'`

    set enddate = `$TIME_STAMP -ehf digital`
    if ( $enddate == "" ) set enddate = tmp`date '+%j%H%M%S'`
    set fyear = `echo $enddate | cut -c -4`

    cd $WORKDIR/rundir
    cat time_stamp.out

########################################################################
# save ascii output files
########################################################################

    if ( ! -d $WORKDIR/ascii ) mkdir $WORKDIR/ascii
    if ( ! -d $WORKDIR/ascii ) then
     echo "ERROR: $WORKDIR/ascii is not a directory."
     exit 1
    endif

    foreach out (`ls *.out *.results input*.nml *_table`)
      mv $out $begindate.$out
    end

    tar cvf - *\.out *\.results | gzip -c > $WORKDIR/ascii/$begindate.ascii_out.tgz

    sbatch --export=source=$WORKDIR/ascii/$begindate.ascii_out.tgz,destination=gfdl:$gfdl_archive/ascii/$begindate.ascii_out.tgz,extension=null,type=ascii --output=$HOME/STDOUT/%x.o%j $SEND_FILE

########################################################################
# move restart files
########################################################################

    cd $WORKDIR

    if ( ! -d $WORKDIR/restart ) mkdir -p $WORKDIR/restart

    if ( ! -d $WORKDIR/restart ) then
      echo "ERROR: $WORKDIR/restart is not a directory."
      exit
    endif

    find $WORKDIR/rundir/RESTART -iname '*.res*' > $WORKDIR/rundir/file.restart.list.txt
    find $WORKDIR/rundir/RESTART -iname '*_data*' >> $WORKDIR/rundir/file.restart.list.txt
    set resfiles     = `wc -l $WORKDIR/rundir/file.restart.list.txt | awk '{print $1}'`

   if ( $resfiles > 0 ) then

      set dateDir = $WORKDIR/restart/$enddate
      set restart_file = $dateDir

      set list = `ls -C1 $WORKDIR/rundir/RESTART`
#      if ( $irun < $segmentsPerJob ) then
#        rm -r $workDir/INPUT/*.res*
#        foreach index ($list)
#          cp $workDir/RESTART/$index $workDir/INPUT/$index
#        end
#      endif

      if ( ! -d $dateDir ) mkdir -p $dateDir

      if ( ! -d $dateDir ) then
        echo "ERROR: $dateDir is not a directory."
        exit
      endif

      foreach index ($list)
        mv $WORKDIR/rundir/RESTART/$index $restart_file/$index
      end

      sbatch --export=source=$WORKDIR/restart/$enddate,destination=gfdl:$gfdl_archive/restart/$enddate,extension=tar,type=restart --output=$HOME/STDOUT/%x.o%j $SEND_FILE

   endif


########################################################################
# move history files
########################################################################

    cd $WORKDIR

    if ( ! -d $WORKDIR/history ) mkdir -p $WORKDIR/history
    if ( ! -d $WORKDIR/history ) then
      echo "ERROR: $WORKDIR/history is not a directory."
      exit 1
    endif

    set dateDir = $WORKDIR/history/$begindate
    if ( ! -d  $dateDir ) mkdir $dateDir
    if ( ! -d  $dateDir ) then
      echo "ERROR: $dateDir is not a directory."
      exit 1
    endif

    find $WORKDIR/rundir -maxdepth 1 -type f -regex '.*.nc'      -exec mv {} $dateDir \;
    find $WORKDIR/rundir -maxdepth 1 -type f -regex '.*.nc.....' -exec mv {} $dateDir \;

    cd $dateDir
      if ( ! -d ${begindate}_nggps3d ) mkdir -p ${begindate}_nggps3d
      mv nggps3d*.nc* ${begindate}_nggps3d
      mv ${begindate}_nggps3d ../.
      if ( ! -d ${begindate}_tracer3d ) mkdir -p ${begindate}_tracer3d
      mv tracer3d*.nc* ${begindate}_tracer3d
      mv ${begindate}_tracer3d ../.
      if ( ! -d ${begindate}_gfs_physics ) mkdir -p ${begindate}_gfs_physics
      mv gfs_physics*.nc* ${begindate}_gfs_physics
      mv ${begindate}_gfs_physics ../.
      if ( ! -d ${begindate}_cloud3d ) mkdir -p ${begindate}_cloud3d
      mv cloud3d*.nc* ${begindate}_cloud3d
      mv ${begindate}_cloud3d ../.

    cd $WORKDIR/rundir

    sbatch --export=source=$WORKDIR/history/$begindate,destination=gfdl:$gfdl_archive/history/$begindate,extension=tar,type=history --output=$HOME/STDOUT/%x.o%j $SEND_FILE
    #sbatch --export=source=$WORKDIR/history/${begindate}_nggps3d,destination=gfdl:$gfdl_archive/history/${begindate}_nggps3d,extension=tar,type=history --output=$HOME/STDOUT/%x.o%j $SEND_FILE
    #sbatch --export=source=$WORKDIR/history/${begindate}_tracer3d,destination=gfdl:$gfdl_archive/history/${begindate}_tracer3d,extension=tar,type=history --output=$HOME/STDOUT/%x.o%j $SEND_FILE
    #sbatch --export=source=$WORKDIR/history/${begindate}_gfs_physics,destination=gfdl:$gfdl_archive/history/${begindate}_gfs_physics,extension=tar,type=history --output=$HOME/STDOUT/%x.o%j $SEND_FILE
    #sbatch --export=source=$WORKDIR/history/${begindate}_cloud3d,destination=gfdl:$gfdl_archive/history/${begindate}_cloud3d,extension=tar,type=history --output=$HOME/STDOUT/%x.o%j $SEND_FILE
