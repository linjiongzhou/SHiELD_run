#!/bin/tcsh -f
#SBATCH --output=/lustre/f2/scratch/Linjiong.Zhou/SHiELD/stdout/%x.o%j
#SBATCH --job-name=C768_20150801.00Z
#SBATCH --partition=batch
#SBATCH --account=gfdl_w
#SBATCH --time=06:00:00
#SBATCH --cluster=c3
#SBATCH --nodes=96
#SBATCH --export=NAME=20150801.00Z,MEMO=_RT2018,EXE=x,LX=16,ALL

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
set CASE = "C768"
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
#set ICS_new = ${INPUT_DATA}/FV3GFS_ICs.v20190701/data/${NAME}_IC
set ICS  = ${INPUT_DATA}/global.v202012/${CASE}/${NAME}_IC
set FIX  = ${INPUT_DATA}/fix.v202104
set GRID = ${INPUT_DATA}/global.v202012/${CASE}/GRID
set FIX_bqx  = ${INPUT_DATA}/climo_data.v201807
set FIX_sfc = ${GRID}/fix_sfc

# sending file to gfdl
set gfdl_archive = /archive/${USER}/SHiELD/${RELEASE}/${NAME}.${CASE}.${TYPE}.${MODE}.${MONO}${MEMO}/
set SEND_FILE =  ${BUILD_AREA}/site/send_file_slurm.csh
set TIME_STAMP = ${BUILD_AREA}/site/time_stamp.csh

# changeable parameters
    # dycore definitions
    set npx = "769"
    set npy = "769"
    set npz = "91"
    set layout_x = $LX
    set layout_y = "16" 
    set io_layout = "1,1"
    set nthreads = "4"

    # blocking factor used for threading and general physics performance
    @ blocksize = ( ${npx} - 1 ) / ${layout_x} * ( ${npy} - 1 ) / ${layout_y}

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
      set n_split = "8"
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
cat ${RUN_AREA}/diag_table_6species_cosp >> diag_table

# copy over the other tables and executable
cp ${RUN_AREA}/data_table data_table
cp ${RUN_AREA}/field_table_6species_aero field_table
cp $executable .

mkdir -p INPUT

# Grid and orography data
ln -sf ${GRID}/* INPUT/

# Date specific ICs
ln -sf ${ICS}/* INPUT/
#ln -sf ${ICS_new}/*sfc_data* INPUT/

# aerosol data
if ( $io_layout == "1,1" ) then
	ln -sf /lustre/f2/dev/gfdl/Linjiong.Zhou/fvGFS_INPUT_DATA/MERRA2/$CASE/*.nc INPUT/
else
	ln -sf /lustre/f2/dev/gfdl/Linjiong.Zhou/fvGFS_INPUT_DATA/MERRA2/$CASE/*.nc.* INPUT/
endif

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
       dnats = 2
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
       do_aerosol = .T.
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
       do_cosp        = .true.
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
       do_sedi_heat = .false.
       vi_max = 1.
       vs_max = 2.
       vg_max = 12.
       vr_max = 12.
       prog_ccn = .true.
       tau_l2v = 225.
       dw_land = 0.16
       dw_ocean = 0.10
       ql_mlt = 1.0e-3
       qi0_crt = 8.0e-5
       rh_inc = 0.30
       rh_inr = 0.30
       rh_ins = 0.30
       ccn_l = 300.
       ccn_o = 200.
       c_paut = 0.5
       c_pracw = 0.8				! aero
       c_psaci = 0.05				! aero
       !c_pracw = 0.35				! aero_cpsd
       !c_psacw = 1.0				! aero_cpsd
       !c_pgacw = 1.e-4				! aero_cpsd
       !c_praci = 1.0				! aero_cpsd
       !c_psaci = 0.35				! aero_cpsd
       !c_pgaci = 0.05				! aero_cpsd
       do_cld_adj = .true.
       use_rhc_revap = .true.
       f_dq_p = 3.0
       rewmax = 10.0
       rermin = 10.0
       !do_new_acc_water = .true.	! aero_cpsd
       !do_psd_water_fall = .true.	! aero_cpsd
       !n0w_sig = 1.2				! aero_cpsd
       !n0w_exp = 66				! aero_cpsd
       !muw = 11.0					! aero_cpsd
       !alinw = 3.e7				! aero_cpsd
       !blinw = 2.0					! aero_cpsd
       !rewflag = 4					! aero_cpsd
       !do_new_acc_ice = .true.		! aero_cpsd
       !do_psd_ice_fall = .true.	! aero_cpsd
       !n0i_sig = 1.0				! aero_cpsd
       !n0i_exp = 10				! aero_cpsd
       !mui = 1.0					! aero_cpsd
       !alini = 11.72				! aero_cpsd
       !blini = 0.41				! aero_cpsd
       !reiflag = 7					! aero_cpsd
       snow_grauple_combine = .false.
/

&COSP_INPUT
  NPOINTS_IT=0,! Max number of gridpoints to be processed in one iteration
  NCOLUMNS=20,  ! Number of subcolumns
  USE_VGRID=.true., ! Use fixed vertical grid for outputs? (if .true. then you need to define number of levels with Nlr)
  NLVGRID=40,       ! Number of levels in statistical outputs (only used if USE_VGRID=.true.)
  CSAT_VGRID=.true., ! CloudSat vertical grid? (if .true. then the CloudSat standard grid is used for the outputs.
                     !  USE_VGRID needs also be .true.)
  !----------------------------------------------------------------------------------
  !--------------- Inputs related to radar simulations
  !----------------------------------------------------------------------------------
  cloudsat_RADAR_FREQ=94.0, ! CloudSat radar frequency (GHz)
  SURFACE_RADAR=0, ! surface=1, spaceborne=0
  cloudsat_use_gas_abs=1,   ! include gaseous absorption? yes=1,no=0
  cloudsat_do_ray=0,        ! calculate/output Rayleigh refl=1, not=0
  cloudsat_k2=-1,           ! |K|^2, -1=use frequency dependent default
  use_precipitation_fluxes=.true.,  ! True if precipitation fluxes are input to the algorithm
  cloudsat_micro_scheme='MMF_v3_single_moment', !'MMF_v3.5_two_moment'
  !----------------------------------------------------------------------------------
  !---------------- Inputs related to lidar simulations
  !----------------------------------------------------------------------------------
  lidar_ice_type=0,    ! Ice particle shape in lidar calculations (0=ice-spheres ; 1=ice-non-spherical)
  OVERLAP=3,           !  overlap assumption used by scops: 1=max, 2=rand, 3=max/rand
  !----------------------------------------------------------------------------------
  !---------------- Inputs related to ISCCP simulator
  !----------------------------------------------------------------------------------
  ISCCP_TOPHEIGHT=1,  !  1 = adjust top height using both a computed
                       !  infrared brightness temperature and the visible
                       !  optical depth to adjust cloud top pressure. Note
                       !  that this calculation is most appropriate to compare
                       !  to ISCCP data during sunlit hours.
                      !  2 = do not adjust top height, that is cloud top
                       !  pressure is the actual cloud top pressure
                       !  in the model
                      !  3 = adjust top height using only the computed
                       !  infrared brightness temperature. Note that this
                       !  calculation is most appropriate to compare to ISCCP
                       !  IR only algortihm (i.e. you can compare to nighttime
                       !  ISCCP data with this option)
  ISCCP_TOPHEIGHT_DIRECTION=2,   ! direction for finding atmosphere pressure level
                                 ! with interpolated temperature equal to the radiance
                                 ! determined cloud-top temperature
                                 ! 1 = find the *lowest* altitude (highest pressure) level
                                 ! with interpolated temperature equal to the radiance
                                 ! determined cloud-top temperature
                                 ! 2 = find the *highest* altitude (lowest pressure) level
                                 ! with interpolated temperature equal to the radiance
                                 ! determined cloud-top temperature. This is the
                                 ! default value since V4.0 of the ISCCP simulator.
                                 ! ONLY APPLICABLE IF top_height EQUALS 1 or 3
  !----------------------------------------------------------------------------------
  !-------------- RTTOV inputs
  !----------------------------------------------------------------------------------
  rttov_Platform=1,    ! satellite platform
  rttov_Satellite=15,  ! satellite
  rttov_Instrument=5,  ! instrument
  rttov_Nchannels=3,   ! Number of channels to be computed
  rttov_Channels=1,2,3,        ! Channel numbers (please be sure that you supply Nchannels)
  rttov_Surfem=0.0,0.0,0.0,  ! Surface emissivity (please be sure that you supply Nchannels)
  rttov_ZenAng=50.0, ! Satellite Zenith Angle
  CO2=5.241e-04, ! Mixing ratios of trace gases
  CH4=9.139e-07,
  N2O=4.665e-07,
  CO=2.098e-07
/

&COSP_OUTPUT
  !- CloudSat
  Lcfaddbze94=.false.,
  Ldbze94=.false.,
  !- CALIPSO
  Latb532=.false.,
  LcfadLidarsr532=.false.,
  Lclcalipso=.false.,
  Lclhcalipso=.true.,
  Lcllcalipso=.true.,
  Lclmcalipso=.true.,
  Lcltcalipso=.true.,
  LparasolRefl=.false.,
  ! CALIPSO phase diagnostics
  Lclcalipsoliq=.false.,
  Lclcalipsoice=.false.,
  Lclcalipsoun=.false.,
  Lclcalipsotmp=.false.,
  Lclcalipsotmpliq=.false.,
  Lclcalipsotmpice=.false.,
  Lclcalipsotmpun=.false.,
  Lclhcalipsoliq=.true.,
  Lcllcalipsoliq=.true.,
  Lclmcalipsoliq=.true.,
  Lcltcalipsoliq=.true.,
  Lclhcalipsoice=.true.,
  Lcllcalipsoice=.true.,
  Lclmcalipsoice=.true.,
  Lcltcalipsoice=.true.,
  Lclhcalipsoun=.true.,
  Lcllcalipsoun=.true.,
  Lclmcalipsoun=.true.,
  Lcltcalipsoun=.true.,
  ! CALIPSO OPAQ diagnostics
  Lclopaquecalipso=.true.,
  Lclthincalipso=.true., 
  Lclzopaquecalipso=.true.,
  Lclcalipsoopaque=.false., 
  Lclcalipsothin=.false.,  
  Lclcalipsozopaque=.false.,
  Lclcalipsoopacity=.false., 
  Lclopaquetemp=.true., 
  Lclthintemp=.true., 
  Lclzopaquetemp=.true., 
  Lclopaquemeanz=.true., 
  Lclthinmeanz=.true., 
  Lclthinemis=.true., 
  Lclopaquemeanzse=.true.,
  Lclthinmeanzse=.true., 
  Lclzopaquecalipsose=.true.,
  ! GROUND LIDAR diagnostics  
  LlidarBetaMol532gr=.false.,  
  LcfadLidarsr532gr=.false.,  
  Latb532gr=.false.,  
  LclgrLidar532=.false.,
  LclhgrLidar532=.false.,  
  LcllgrLidar532=.false., 
  LclmgrLidar532=.false.,
  LcltgrLidar532=.false.,
  ! ATLID diagnostics
  LlidarBetaMol355=.false.,
  LcfadLidarsr355=.false., 
  Latb355=.false., 
  Lclatlid=.false., 
  Lclhatlid=.false., 
  Lcllatlid=.false.,
  Lclmatlid=.false.,
  Lcltatlid=.false.,
  !- ISCCP
  Lalbisccp=.false.,
  Lboxptopisccp=.false.,
  Lboxtauisccp=.false.,
  Lpctisccp=.false.,
  Lclisccp=.false.,
  Ltauisccp=.false.,
  Lcltisccp=.false.,
  Lmeantbisccp=.false.,
  Lmeantbclrisccp=.false.,
  !- MISR
  LclMISR=.false.,
  !- Use lidar and radar
  Lclcalipso2=.false.,
  Lcltlidarradar=.false.,
  Lcloudsat_tcc=.false.,
  Lcloudsat_tcc2=.false.,
  !- These are provided for debugging or special purposes
  Lfracout=.false.,
  LlidarBetaMol532=.false.,  
  !- MODIS
  Lcltmodis=.true.,
  Lclwmodis=.true.,
  Lclimodis=.true.,
  Lclhmodis=.true.,
  Lclmmodis=.true.,
  Lcllmodis=.true.,
  Ltautmodis=.true.,
  Ltauwmodis=.true.,
  Ltauimodis=.true.,
  Ltautlogmodis=.true.,
  Ltauwlogmodis=.true.,
  Ltauilogmodis=.true.,
  Lreffclwmodis=.true.,
  Lreffclimodis=.true.,
  Lpctmodis=.true.,
  Llwpmodis=.true.,
  Liwpmodis=.true.,
  Lclmodis=.false.,
  !- RTTOV
  Ltbrttov=.false.,
  ! -CLOUDSAT precipitation frequency/occurence diagnostics
  Lptradarflag0=.false.,
  Lptradarflag1=.false.,
  Lptradarflag2=.false.,
  Lptradarflag3=.false.,
  Lptradarflag4=.false.,
  Lptradarflag5=.false.,
  Lptradarflag6=.false.,
  Lptradarflag7=.false.,
  Lptradarflag8=.false.,
  Lptradarflag9=.false.,
  Lradarpia=.false.,
  !- CloudSat+MODIS joint diagnostics
  Lwr_occfreq=.false.,
  Lcfodd=.false.
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
