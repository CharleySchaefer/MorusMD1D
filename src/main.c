#include "main.h"

double get_center_of_mass(double *R, int n_beads){
  int i;
  double Rcm; /* centre of mass */

  if(n_beads>2) {
    Rcm=0.0; 
    for (i=1; i<n_beads-1; i++)
      Rcm+=R[i];
    Rcm/=(n_beads-2);  
  
  } else{        /* Dumbbell model */
    Rcm=0.5*(R[n_beads-1]+R[0]);
  }
  return Rcm;
}

/*
  Force in          [kT * unit length]
  Eact in units     [kT]
  diss_length in    [unit length]

  delta_Eact = min(   Eact_diss  ,   |force|*lengthscale  )
*/
double get_sticker_diss_coeff(double force, double Eact_diss, double sticker_diss_length){

  /* |force|*lengthscale */
  double delta_Eact=(force>0?force:-force)*sticker_diss_length;

  /* minimum */
  delta_Eact=(delta_Eact<Eact_diss?delta_Eact:Eact_diss);
  return exp(delta_Eact);
}



int main(int argc, char *argv[])
{
  /* Declare arrays */
  double *hist_mu=NULL, *hist_mu_all=NULL, *hist_mu_all_count=NULL, *R=NULL, *dR=NULL, *dRdi=NULL, *dRdt=NULL, *ks=NULL, *lambda_loc;
    /* Lihktman's method to calculate mu(t) */
  int    Nx=-1;                  /* discretise space: number of grid points*/
  int    xiL, xiU;               /* lattice indices of chain ends          */
  int    xiL0, xiU0;             /*  "         "    in previous time step  */
  double lambda_max=-1;          /* finite-extensibility only if lambda_max>0  */
  double Lx,                     /* discretise space: interval length      */
         dx,                     /* discretise space: lattice spacing      */ 
         *mu_instant=NULL, *mu_LMcumm=NULL,/* mu_instant(t) and cummulative mu(t)    */ 
         *t_alive=NULL;               /* 'lifetime' of grid segment             */
  int    Ncumm=0;                /* number of counts in cummulative mu_LM  */
  int    *hist_count=NULL, *state=NULL,
         *sequence=NULL;         /* read using --sequence <string>
                                    length:  n_beads 
                                    value 0: bead is not a sticker
                                    value 1: bead is a sticker */ 
  char *arg=NULL; 

  /*===========================================*/
  /* DEFAULT USER SETTINGS                     */
  int new_simulation=1;          /* 1: create new output file; 0(TODO): read previous output */
  char opt_boundary='B';
  int cheat_shear=0;

  double lambda0=1.0; /* Initial (t=0) stretch ratio */

  /* Chain length */
  int      Ze=-1;                /* number of entanglements */
  int      n_beads=-1;           /* number of beads */
             
  /* Time */
  double   ttime=0.0,           /* time */ 
           ttopen=0.0,           /* total time sticker 1 is open - to validate correct value of p */  
           time_fac=0.4;         /* step size; must be >0 and <0.5 */

  /* Flow */
  double   Weissenberg=0.0;
  long int Ntime = 2000000, 
           Nprint=    1000;      /* overriden by --export-logt */
 
  /* Thermal fluctuations  */
  int rseed=0;                   /* Random number seed. See --rseed argument */
  int include_fluctuations=1;

    
  /* Stickers */
  int      Zs=0;                 /* number of stickers per chain */
  int sync_stickers=0;           /* 1: stickers are all open or all closed */

  /* Parameters for strain-enhanced sticker dissociation */
  int accelerated_sticker_dissociation=0;
  double Eact_diss=0;             /* 0: no influence of tension in the chain         */
  double sticker_diss_length=0.1; /* in nm                                           */ 
  double b_Kuhn=0.36;             /* in nm  (0.36 for amino acid)                    */
  double Ne=552.5;                /* length of entanglement strand [# of Kuhn lengths] (552.5 for silk)*/
  double force_fac;               /* calculated during initialisation*/
  double sticker_diss_coeff;      /* calculated at every time step                   */
  double sticker_stress;

  double sticker_kbondswap=0.0,  /* sticker bondswap rate */
         sticker_kdiss=0.0,      /* sticker dissociation rate */
         sticker_p=0.0;          /* mean fraction of closed stickers */  

  /* Output */
  int opt_export_sequence=0;
  int opt_export_logt=0; /*0: export every Nprint time steps; 
                           1: use logarithmic interval */

     /* Output: mu(t) */
  int opt_calc_mu=0,
      mu_Nanalyse=1, /* update bin values of mu(t) every mu_Nanalyse time steps */
      mu_nbin=200;   /* number of bins for time interval */
  double mu_tL=1e-5, /* lower limit     of time interval */
         mu_tU=1e7;  /* upper limit     of time interval */

  /* END DEFAULT USER SETTINGS                 */
  /*===========================================*/


  /*===========================================*/
  /* PARSE PROGRAM ARGUMENTS*/
  int i, j, check;
  char *outdir=(char*)malloc(1000*sizeof(char));
  outdir[0]='.'; outdir[1]='\0'; /* default */
  char *file_timeprogress=(char*)malloc(1000*sizeof(char));
  char *file_param=(char*)malloc(1000*sizeof(char));
  FILE *ifp;
  
  if(argc==1){
    printf("No input given. For input options, run %s --help\n", argv[0]); goto  TERMINATE_PROGRAM;
  }  
  i=1; while( i < argc ){

    arg = argv[i];
    if ( strcmp( arg, "--help" ) == 0 ){ 
      printHelp(argv[0]);
      goto  TERMINATE_PROGRAM;
    } else if ( strcmp( arg, "--version" ) == 0 ){ /* Random number seed */
      printf("%s %s\n", argv[0]+2, VERSION);  exit(1);
      i += 1;
    } else if ( strcmp( arg, "--outdir" ) == 0 ){ /* dir */
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        sprintf(outdir, "%s", argv[i+1]); i+=2; 
        printf("  Output directory:          %s\n", outdir);
      }
    } else if ( strcmp( arg, "--time-start" ) == 0 ){ /* Random number seed */
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        ttime=atof(argv[i+1]); i+=2;
        printf("  Time at start:                %e\n", ttime);
      }
    } else if ( strcmp( arg, "--export-logt" ) == 0 ){ 
      opt_export_logt=1;
      i+=1;
/* CALCULATION OF mu(t)  */
    } else if ( strcmp( arg, "--export-mu" ) == 0 ){ 
      opt_calc_mu=1;
      i+=1;
    } else if ( strcmp( arg, "--mu-time-steps" ) == 0 ){
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        mu_Nanalyse=atoi(argv[i+1]); i+=2;
        opt_calc_mu=1; 
        printf("  Analyse mu(t) after every %d time steps.\n", mu_Nanalyse);
      }
    } else if ( strcmp( arg, "--mu-Nbins" ) == 0 ){ 
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        mu_nbin=atoi(argv[i+1]); i+=2;
        opt_calc_mu=1; 
        printf("  Number of bins for mu(t) histogram: %d\n", mu_nbin);
      }
    } else if ( strcmp( arg, "--mu-time-L" ) == 0 ){ 
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        mu_tL=atof(argv[i+1]); i+=2;
        opt_calc_mu=1; 
        printf("  Lower time limit for mu(t) histogram: %e\n", mu_tL);
      }
    } else if ( strcmp( arg, "--mu-time-U" ) == 0 ){
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        mu_tU=atof(argv[i+1]); i+=2;
        opt_calc_mu=1; 
        printf("  Upper time limit for mu(t) histogram: %e\n", mu_tU);
      }
    } else if ( strcmp( arg, "--mu-Nx" ) == 0 ){
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        Nx=atoi(argv[i+1]); opt_calc_mu=1; i+=2;
        printf("  Discretisation of space - Nx: %d\n", Nx);
      }
    } else if ( strcmp( arg, "--mu-dx" ) == 0 ){
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        dx=atof(argv[i+1]); opt_calc_mu=1; i+=2;
        printf("  Discretisation of space - dx: %f\n", dx);
      }
/* THERMAL FLUCTUATIONS */
    } else if ( ( strcmp( arg, "--exclude-fluctuations" ) == 0 ) ||
                ( strcmp( arg, "--unset-fluctuations"   ) == 0 )  ){ 
      include_fluctuations=0;
      i+=1;
    } else if ( strcmp( arg, "--export-sequence" ) == 0 ){ 
      opt_export_sequence=1;
      i+=1;
    } else if ( strcmp( arg, "--set-boundary" ) == 0 ){ /* Boundary - A,B or C*/
      if( i+1>=argc ) {
         printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        opt_boundary=argv[i+1][0]; i+=2;
     }
    } else if ( strcmp( arg, "--init-continue" ) == 0 ){ 
      new_simulation=0;
      printf("Error: --init-continue not yet implemented.\n"); 
      i+=1;
    } else if ( strcmp( arg, "--Wi" ) == 0 ){ /* Weissenberg number */
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        Weissenberg=atof(argv[i+1]); i+=2;
        printf("  Weissenberg number:           %f\n", Weissenberg);
      }
    } else if ( strcmp( arg, "--cheat-shear" ) == 0 ){ 
      cheat_shear=1;
      i+=1;
    } else if ( strcmp( arg, "--Ntime" ) == 0 ){ /* Random number seed */
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        Ntime=(long int)atol(argv[i+1]); i+=2;
        printf("  Number of time steps:         %ld\n", Ntime);
      }
    } else if ( strcmp( arg, "--time-fac" ) == 0 ){ /* Random number seed */
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        time_fac=atof(argv[i+1]); i+=2;
        if(time_fac>=0.5 || time_fac<=0) {
          printf("Error: --time-fac <fac> should have a positive value smaller than 0.5.\n"); goto  TERMINATE_PROGRAM;
        }
        printf("  Time step prefactor:          %f\n", time_fac);
      }
    } else if ( strcmp( arg, "--Nprint" ) == 0 ){ /* Random number seed */
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        Nprint=(long int)atol(argv[i+1]); i+=2;
        printf("  Print after time steps:       %ld\n", Nprint);
      }
    } else if ( strcmp( arg, "--rseed" ) == 0 ){ /* Random number seed */
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        rseed=atoi(argv[i+1]); i+=2;
        printf("  Random number seed:           %d\n", rseed);
      }
    } else if ( strcmp( arg, "--Ze" ) == 0 ){ /* Number of entanglements*/
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        Ze=atoi(argv[i+1]); i+=2;
        printf("  Number of entanglements (Ze): %d\n", Ze);
      }
    } else if ( strcmp( arg, "--Ne" ) == 0 ){ /* Number of entanglements*/
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        Ne=atof(argv[i+1]); i+=2;
        printf("  Monomers per entanglement strand (Ne): %f\n", Ne);
      }
    } else if ( strcmp( arg, "--kuhn-length" ) == 0 ){ /* Number of entanglements*/
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        b_Kuhn=atof(argv[i+1]); i+=2;
        printf("  Kuhn length: %f\n", b_Kuhn);
      }
    } else if ( strcmp( arg, "--lambda-max" ) == 0 ){ /* Number of entanglements*/
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        lambda_max=atof(argv[i+1]); i+=2;
        printf("  Maximum stretch ratio lambda_max: %e\n", lambda_max);
      }
    } else if ( strcmp( arg, "--sticker-sync") == 0 ){ /* synchronise stickers: either all open or all closed*/
        sync_stickers=1; i++;
    } else if ( strcmp( arg, "--sticker-p" ) == 0 ){ /* */
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        sticker_p=atof(argv[i+1]); i+=2;
        printf("  Fraction of closed stickers:        %e\n", sticker_p);
      }
    } else if ( strcmp( arg, "--sticker-kbondswap" ) == 0 ){ /* */
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        sticker_kbondswap=atof(argv[i+1]); i+=2;
        printf("  Sticker bondswap rate:        %e\n", sticker_kbondswap);
      }
    } else if ( strcmp( arg, "--sticker-kdiss" ) == 0 ){ /* */
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;  
      } else {
        sticker_kdiss=atof(argv[i+1]); i+=2;
        printf("  Sticker dissociation rate:        %e\n", sticker_kdiss);
      }
    } else if ( strcmp( arg, "--sticker-kdiss-length" ) == 0 ){ /* */
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;  
      } else {
        sticker_diss_length=atof(argv[i+1]); i+=2;
        printf("  Sticker dissociation length:        %e\n", sticker_diss_length);
      }
    } else if ( strcmp( arg, "--sticker-kdiss-Eact" ) == 0 ){ /* */
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;  
      } else {
        Eact_diss=atof(argv[i+1]); i+=2;
        printf("  Activation energy of sticker dissociation:        %e\n", Eact_diss);
      }
    } else if ( strcmp( arg, "--sequence" ) == 0 ){ /* Random number seed */
      if (i+1>=argc) {
        printf("Error: unexpected program arguments.\n"); goto  TERMINATE_PROGRAM;
      } else {
        n_beads=strlen(argv[i+1]);
        sequence=(int*)calloc(n_beads, sizeof(int));
        for(j=0; j<n_beads; j++) {
          if( '1'==argv[i+1][j] ){ /* if bead j is a sticker */
            Zs         +=1;
            sequence[j] =1;
          }
        }
        printf("  Number of beads:              %d\n", n_beads);
        printf("  Number of stickers      (Zs): %d\n", Zs);
        i+=2;
      }
    } else {
      printf( "Error: Argument \"%s\" not recognised!\n", arg );
      printf( "       For help, run \"%s --help\".\n", argv[0] );
      
      goto  TERMINATE_PROGRAM;
    }
  }
  srand(rseed);


  /* */
  if(opt_export_logt==1) {
    Nprint=1;
  }

 

  /* use file_timeprogress as buffer to create directory*/
  sprintf(file_timeprogress, "mkdir -p %s",      outdir);
  system( file_timeprogress); 
  /* define file_timeprogress */
  sprintf(file_timeprogress, "%s/%s", outdir, "timeprogress.out"); 
printf("output: %s\n", file_timeprogress);
  sprintf(file_param,        "%s/%s", outdir, "param.in");
  ifp=fopen(file_param, "w");
  i=0; while( i < argc ){
    fprintf(ifp, "%s ", argv[i]);
    i++;
  }
  fprintf(ifp, "\n");
  fclose(ifp);


  if(Ze==-1) {
    printf("Error: argument --Ze missing.\n."); goto  TERMINATE_PROGRAM;
  }



  /* END PARSING PROGRAM ARGUMENTS             */
  /*===========================================*/


  /*===========================================*/
  /* DECLARATIONS  */
  long int    ti;           /* time step index         */
  double logt;              /* logarithm of time       */
  double mu, maxR0, minRZ;  /* used to calculate mu(t) */

  double mu_logtL=log10(mu_tL),
         mu_logtU=log10(mu_tU),
         mu_dlogt=(mu_logtU-mu_logtL)/(mu_nbin-1),
         hist_binwidth=mu_dlogt; /* TODO: replace all hist_binwidth instances by mu_dlogt  */


  double ks_norm; /* for finite chain extensibility */
  //double lambda_loc;
  if(lambda_max>0){
     ks_norm=get_ks_norm(lambda_max);
  }

  if(Eact_diss>0) {
    accelerated_sticker_dissociation=1;
    printf("  Stretch-enhanced sticker dissociation activated with parameters:\n");
    printf("    Ne:                  %f\n", Ne);
    printf("    Eact:                %f\n", Eact_diss);
    printf("    Dissociation length: %f\n", sticker_diss_length);
    printf("    Kuhn length:         %f\n", b_Kuhn);
  }

  double strain_rate=Weissenberg/(Ze*Ze);
  printf("  Strain rate: %f\n", strain_rate);

  double Rcm=0.0, /* centre of mass */ 
         stress;  /* chain stress   */
  /*===========================================*/


  /*===========================================*/
  /* MEMORY ALLOCATION  */
  if(Nx>0) {
    printf("allocating memory for mu(t) calculation.\n");
    t_alive     = (double*)calloc(Nx, sizeof(double));
    mu_instant  = (double*)calloc(mu_nbin, sizeof(double));
    mu_LMcumm   = (double*)calloc(mu_nbin, sizeof(double));
  }
  hist_mu     = (double*)calloc(mu_nbin,sizeof(double));
  hist_mu_all = (double*)calloc(mu_nbin,sizeof(double));
  hist_mu_all_count = (double*)calloc(mu_nbin,sizeof(double));
  R     = (double*)calloc(n_beads, sizeof(double)); /* coordinate */
  dR  = (double*)calloc(n_beads, sizeof(double)); /* local stretch  */
  dRdi  = (double*)calloc(n_beads, sizeof(double)); 
  lambda_loc=(double*)calloc(n_beads, sizeof(double));/* local stretch ratio */
  ks  = (double*)calloc(n_beads, sizeof(double)); /* non-linearity coefficient of spring constant */
  dRdt  = (double*)malloc(n_beads* sizeof(double));
  hist_count  = (int*)calloc(mu_nbin,sizeof(int));  /* mu(t) histogram */
  state = (int*   )calloc(n_beads, sizeof(int));    /* state of a bead */
  /*===========================================*/


  /*===========================================*/
  /* INITIALISE                                */
  int count_relaxations=0;

  /* discretisation */
    /* spatial discretisation */
  double di;
  /* 
     Boundary condition A:
     Real chain:          E-------------------E   ; E=chain end
     Bead distribution:   B---B---B---B---B---B   ; B=bead
     Postion:             a   b   c 
                        a= -Ze/2
                        b= -Ze/2 + di
                        c= -Ze/2 + 2*di, etc.

     NOTE: In this implementation the distance a-b is fixed; this 
           is slightly different from that used by Alexei Lihktman, who
           used 
                              E-------------------E      ; E=chain end
                          C---B---B---B---B---B---B---C  ; B=bead; C=implicit bead

           where C the distance C-B is fixed, and C is not explicitly
           counted as a bead. 

     Boundary condition B:
     Real chain:          E-------------------E   ; E=chain end
     Bead distribution: B---B---B---B---B---B---B ; B=bead
     Postion:           a b c   d  
                        a= -Ze/2 - dcorr ; dcorr=di/2
                        b= -Ze/2
                        c= -Ze/2 + di
                        d= -Ze/2 + 2*di, etc.
  */
  switch (opt_boundary){
    case 'A': { di = (double)Ze/(n_beads-1);                         break;}
    case 'B': { di = (double)Ze*(n_beads-1)/(n_beads-2)/(n_beads-1); break;}
    case 'C': { di = (double)Ze/(n_beads-1);                         break;}
    default:  { printf("Error: unknown boundary option %c.\n", opt_boundary);
                goto TERMINATE_PROGRAM;}
  }
  

  double id2i=1.0/( M_PI*M_PI*di*di ); 
  for (i=0; i<n_beads; i++){
    R[i] = (i-(n_beads-1)*0.5)*di;    /* spatial coordinate of each bead                  */
    ks[i] = 1.0;                      /* non-linearity coefficient of the spring constant */
  }
  if(accelerated_sticker_dissociation){
    force_fac=3*id2i*M_PI*M_PI/(b_Kuhn*pow(Ne, 1.5)); 
    printf("numerical force_fac: %e\n", force_fac);
  }

  double   dcorr=( R[n_beads-1]-R[0]-Ze )*0.5 ;
  maxR0=R[0]+dcorr;         /* chain end: Rmin=R[0        ]+dcorr */
  minRZ=R[n_beads-1]-dcorr; /* chain end: Rmax=R[n_beads-1]-dcorr */
  double TubeLength0=minRZ-maxR0;
  double time_L0=0.0;
    /* time discretisation */
  double dt=time_fac/id2i; /* */
// for (i=0; i<n_beads; i++)
//    R[i] = 10*(i-(n_beads-1)*0.5)*di;

  /*-----------------------------------------*/
  /* STICKER PROPERTIES */
  double kS_open  = sticker_kdiss + sticker_kbondswap*(1-sticker_p);
  if (sticker_p==1.0){
    kS_open=0.0; 
    printf("Warning: for p=1 the sticker opening rate is set to zero.\n");
   }
  double kS_close = sticker_kdiss*sticker_p/(1-sticker_p) + sticker_kbondswap*sticker_p; 
  kS_open  *= dt; /* on rate  in units of time step */
  kS_close *= dt; /* off rate in units of time step*/

  int countOpen=0,countClose=0;
  double tmpf, Rtmp1, Rtmp2;
  if (Zs>0){
    printf("  Mean fraction of closed stickers: %e\n", sticker_p);
    if(sticker_p==0) /* all stickers are open */
      kS_close=0.0;
    else {
      printf("  Sticker closing probability / time step: %e\n", kS_close);
      printf("  Sticker opening probability / time step: %e\n", kS_open);
      if(Zs>0 && kS_close>=1.0 && sticker_p<1.0){
        printf("Error: probability of sticker closing >1 per time step.\n");
        goto  TERMINATE_PROGRAM;
      }
    }
  }

  /* Close stickers with probability 'sticker_p' */
  if(sticker_p>0) {
    for(i=0; i<n_beads; i++){                             /* loop all mmonomers */
      if(sequence[i]==1)                                  /* if monomer is sticker*/
        if( (double) rand()/(RAND_MAX+1.0) < sticker_p )  /* close sticker with probability */
          state[i] = 1;                                   /* 0: open sticker; 1: closed sticker */
    }
  }



    /* noise amplitude */
  double sigma=sqrt( 2.0/(3.0*M_PI*M_PI*di*dt) ); /* std of noise */
  double uniform_width = sigma*sqrt(12); /* width of uniform distribution
                                            to generate correct std   */

  /* For mu(t) calculation: initialise periodic box and get grid indices  */
  mu=1.0; /* mu at t=0 */
  if(Nx>0){
    /* discretisation of space: periodic at x=0 and at x=L */
    Lx=Nx*dx;
    if(Lx<2*Ze){
      printf("Error: grid length %f should be larger than twice the quiescent chain length %d.\n", Lx, Ze); goto TERMINATE_PROGRAM;
    }
    printf("  Grid length: %f\n", Lx);

    /* Grid indices of chain ends */ 
    xiL0=get_grid_index( R[0        ]+dcorr, Nx, dx );
    xiU0=get_grid_index( R[n_beads-1]-dcorr, Nx, dx ); 
  }


  /* Create output file (header and initial values */
  if(new_simulation) { 
    printf("Create %s\n", file_timeprogress);
    ifp=fopen(file_timeprogress, "w");
    /*header line*/
    fprintf(ifp, "%12s %12s %12s %12s %12s", "time/tauE", "Rmin/a", "Rmax/a", "Rcm/a", "stress");
    if(opt_export_sequence){
      for(i=0; i<n_beads; i++)
        fprintf(ifp, "   R%02d_%04d/a", sequence[i], i+1); /* R<monomer type>_<monomer index>/a */
    }
    fprintf(ifp, "\n"); /* end of header line*/

    Rcm=0.0;     /* centre of mass */
    stress=0.0;  /* calculate stress */
    for (i=0; i<n_beads; i++) {
      Rcm+=R[i];
      if(i>0)
        stress+=pow(R[i]-R[i-1], 2)/di; /* stress */
    }
    Rcm/=(n_beads); 

    /*data line*/
    fprintf(ifp, "%12.5e %12.5e %12.5e %12.5e %12.5e", ttime, R[0]+dcorr-Rcm, R[n_beads-1]-dcorr-Rcm, Rcm, stress);

    if(opt_export_sequence)
      for(i=0; i<n_beads; i++)
        fprintf(ifp, " %12.5e", R[i]-Rcm);
    fprintf(ifp, "\n"); /* end of data line*/

    fclose(ifp);
  }

  

  int n_first=1, n_last=n_beads-2; /* beads that move independently */
  if (opt_boundary=='C') {         
    n_first=0;
    n_last=n_beads-1;
  }  

/*================================================================*/
/* TIME LOOP */
int is_sticker, is_closed;
for(ti=0; ti<Ntime; ti++) {


  /* Get centre of mass (needed for extensional flow) */
  if(strain_rate>0) {
    Rcm=get_center_of_mass(R, n_beads);
  }

  /* Get local stretch ; the stretch ratio at i is given by dR[i]/di */
  Rtmp1=R[0];
  for (i=1; i<n_beads; i++) {
    Rtmp2=R[i];
    dR[i-1]=Rtmp2-Rtmp1;  /* Distance between beads  (dR[i]/di is the local stretch ratio lambda)             */
    if(lambda_max>0) {    /* if finite extensibility is taken into account */
      lambda_loc[i-1]=dR[i-1]/di;
      ks[i-1]=get_ks(lambda_loc[i-1], lambda_max, ks_norm);  
//if(lambda_loc>10)
//  printf("test: %f %f\n", lambda_loc, ks[i-1]);
    }
    Rtmp1=Rtmp2;
  }
  

  /* Get Rate of change  */
  for (i=n_first; i<=n_last; i++) { // only independently moving beads (chain end positions may be dependent)

    /* determine if bead is a sticker, and if it is closed*/
    is_sticker=(sequence[i]==1);
    if(is_sticker) {
      is_closed=(state[i]==1); 




      if(accelerated_sticker_dissociation) { /* calculate opening rate kS_open */

        if (opt_boundary=='C'){    
          if(i==0) {
            sticker_stress = force_fac*( -di + dR[0]  )*ks[0];
          } else if (i==n_beads-1) {
            sticker_stress = force_fac*( -dR[n_beads-2] + di )*ks[n_beads-2];
          } else {
            sticker_stress = force_fac*( ks[i]*dR[i]-ks[i-1]*dR[i-1]  );
          }    
          sticker_diss_coeff=get_sticker_diss_coeff(sticker_stress, Eact_diss, sticker_diss_length);

          kS_open  = sticker_kdiss*sticker_diss_coeff + sticker_kbondswap*(1-sticker_p);
          kS_open *= dt; /* off rate in units of time step*/


        } else {
          if( lambda_loc[i-1]>lambda_max || lambda_loc[i]>lambda_max ) { /* open sticker if lambda becomes too large */
            kS_open=2;
          } else {
            sticker_stress = force_fac*( ks[i]*dR[i]-ks[i-1]*dR[i-1]  );

            sticker_diss_coeff=get_sticker_diss_coeff(sticker_stress, Eact_diss, sticker_diss_length);

            kS_open  = sticker_kdiss*sticker_diss_coeff + sticker_kbondswap*(1-sticker_p);
            kS_open *= dt; /* off rate in units of time step*/

          }
        }
      }
    }

    /* Retraction */
    if(is_sticker && is_closed){ /* closed sticker cannot retract*/
      dRdt[i]=0.0; 
    } else {                      /*  bead free to retract */
      if (opt_boundary=='C'){    
        if(i==0) {
          dRdt[i] = id2i*( -di + dR[0]  )*ks[0];
        } else if (i==n_beads-1) {
          dRdt[i] = id2i*( -dR[n_beads-2] + di )*ks[n_beads-2];
        } else {
          dRdt[i] = id2i*( ks[i]*dR[i]-ks[i-1]*dR[i-1]  );
        }    
      } else {
        dRdt[i] = id2i*( ks[i]*dR[i]-ks[i-1]*dR[i-1]  );
      }
    }

    /* Extension rate */
    if(strain_rate>0) {
      if (cheat_shear==1) { /* 'cheat' shear flow */
        if(             (R[i]-Rcm)<0.5 ) {
          dRdt[i] -= strain_rate;             /* rate constant for distances beyond tube diameter  */
        } else if(      (R[i]-Rcm)>0.5 ) {
          dRdt[i] += strain_rate;             /* rate constant for distances beyond tube diameter  */
        } else {
          dRdt[i] += (R[i]-Rcm)*strain_rate;  /* rate 'extensional' within the tube diameter       */
        }
      } else {              /* extensional flow */
        dRdt[i] += (R[i]-Rcm)*strain_rate;
      }
    }

    /* Thermal fluctuations*/
    if(include_fluctuations && (!is_sticker || !is_closed)) {
      dRdt[i] += uniform_width*(  (double)rand()/(RAND_MAX+1.0) - 0.5   );
    }

    /* Change rate of the sticker Sticker interactions */
    if(!sync_stickers && is_sticker) { /* 0: no sticker; 1: sticker */
      /* sticker is open */
      if(is_closed ) { /* change state: open sticker */
        if ( (double) rand()/(RAND_MAX+1.0) < kS_open  ) {
          countOpen++;
          state[i]=0;
        }
      } else {         /* change state: close sticker */
        if ( (double) rand()/(RAND_MAX+1.0) < kS_close  ) {
          state[i]=1;
          countClose++;
        }
      }
    }
  } /* Loop over monomers done (dRidt calculated) */

  if( sync_stickers ){ /*synced stickers all opened or closed simultaneously */
    if(is_closed){ /* stickers are currently closed */
        if ( (double) rand()/(RAND_MAX+1.0) < kS_open  ) {
          countOpen++;
          if(n_beads>2){
            for (i=n_first; i<=n_last; i++) {
              if( sequence[i]==1)
                state[i]=0;   
            }
          } else {
            state[0]=0;state[n_beads-1]=0; 
          }
        }
    } else {         /* change state: close sticker */
        if ( (double) rand()/(RAND_MAX+1.0) < kS_close  ) {
          countClose++;
          if(n_beads>2){
            for (i=n_first; i<=n_last; i++) {
              if( sequence[i]==1)
                state[i]=1;
            }
          } else {
            state[0]=1;state[n_beads-1]=1; 
          }
        }
    }
  } /* stickers synchronised */


  /* Update time and configuration */
  if(state[0]==0)
    ttopen+=dt;
  ttime+=dt;
  if( opt_boundary=='C' ) {
    for (i=0; i<n_beads; i++) {
      R[i] +=dt*( dRdt[i]   ); 
    }
  } else { /* boundary A or B */
    if(n_beads>2) {
      for (i=1; i<n_beads-1; i++) {
        R[i] += dt*( dRdt[i]   ); 
      }
      R[0]         = R[1]        -di; /* boundary condition at i=0         */
      R[n_beads-1] = R[n_beads-2]+di; /* boundary condition at i=n_beads-1 */
    } else { //n_beads=2  DUMBBELL MODEL

        /* determine if bead is a sticker, and if it is closed*/
      is_sticker=(sequence[i]==1);
      if(is_sticker)
        is_closed=(state[i]==1); 

      /* initialise */
      dRdt[0        ]=0; 
      dRdt[n_beads-1]=0;


      /* extension */
      if(strain_rate>0){
        dRdt[0]         += (R[0        ]-Rcm)*strain_rate;
        dRdt[n_beads-1] += (R[n_beads-1]-Rcm)*strain_rate;
      }
        /* thermal fluctuations */
      if( (!is_sticker || !is_closed)) {
        /* retraction */
        dRdt[0        ] += (0.5/(di*di))*( R[n_beads-1]-R[0] - Ze );
        dRdt[n_beads-1] += -dRdt[0        ];
        /* thermal fluctuations */
        if(include_fluctuations) {
          dRdt[0        ] += sqrt( 0.5*12*2.0/(3.0*di*dt) )*(  (double)rand()/(RAND_MAX+1.0) - 0.5   );
          dRdt[n_beads-1] += sqrt( 0.5*12*2.0/(3.0*di*dt) )*(  (double)rand()/(RAND_MAX+1.0) - 0.5   );
        }
      }

         /* Update */
      R[0        ] += dt*( dRdt[0        ]   ); 
      R[n_beads-1] += dt*( dRdt[n_beads-1]   ); 
      
    }
  }

  /*================================================================*/
  if(opt_calc_mu) {  
    if ((ti%mu_Nanalyse==0) || (ti==Ntime-1)){         /* analyse mu(t) */

      /* METHOD 1 */
      if(Nx>0) {
      xiL=get_grid_index( R[0        ]+dcorr, Nx, dx );
      xiU=get_grid_index( R[n_beads-1]-dcorr, Nx, dx );

      /* Reset lifetime of grid segments in events of end-group retraction*/
      if( xiL<xiL0-Ze/dx) {  /* bead retracted through periodic boundary */
        for( i=xiL0; i!=xiL; i) {
          t_alive[i]=0.0; i++; i=(i==Nx?0:i); /* reset lifetime grid segments*/
        }
      } else if( (xiL>xiL0) && (xiL<xiL0+Ze/dx)) {
        for( i=xiL0; i<xiL; i++) 
          t_alive[i]=0.0;
      }

      if( xiU>xiU0+Ze/dx) {  /* bead retracted through periodic boundary */
        for( i=xiU0; i!=xiU; i) {
          t_alive[i]=0;i--; i=(i==-1?Nx-1:i); /* reset lifetime grid segments*/
        }
      } else if( (xiU<xiU0) && (xiU>xiU0-Ze/dx)) {
        for( i=xiU0; i>xiU; i--) 
          t_alive[i]=0.0;
      }
      /* get mu_instant */
      j=get_mu_instant(mu_instant, mu_nbin, mu_logtL, mu_dlogt, t_alive, Nx, xiL, xiU);    
       /*j = time index corresponing to longest life time of a grid segment*/
      /* save */
      Ncumm++;
      for(i=0; i<j; i++  ) {
        mu_LMcumm[i]+=mu_instant[i];
      }
      /* update lifetimes of alive grid segments */
      xiL0=xiL; xiU0=xiU;
      for( i=xiL-1 ; i!=xiU ; i ){
        i++; i=(i==Nx?0:i); 
        t_alive[i]+=dt;
      }
      }

     

      /* METHOD 2 */
      maxR0=(maxR0>R[0]        +dcorr ? maxR0 : R[0]        +dcorr);
      minRZ=(minRZ<R[n_beads-1]-dcorr ? minRZ : R[n_beads-1]-dcorr);
      mu=fmin(mu, (minRZ-maxR0)/TubeLength0); /* mu(t): fraction of tube unrelaxed since t0 */
      logt=log10(ttime-time_L0);              /* log(t-t0) since start of mu( t-t0 ) measurement  */

      i=(int)((logt-mu_logtL)/hist_binwidth); /* bin number */
      j=i;

      if(mu<=0 || (ti==Ntime-1)) { /* RESET */
      //  for(j=i; j<mu_nbin; j++)
      //    hist_mu_all_count[j]++;
        count_relaxations++;
        time_L0=ttime;             /* time at which mu(t) measurement starts */
        mu=1.0;

        maxR0=R[0]+dcorr; 
        minRZ=R[n_beads-1]-dcorr;
        TubeLength0=minRZ-maxR0;
      
        for(i=0; i<mu_nbin; i++) { /* reset histogram */
          if(/*hist_mu[i]>0 && */ hist_count[i]>0) {
            hist_mu_all[i]      +=hist_mu[i]/hist_count[i];
            hist_mu_all_count[i]++;
          }
          else if (i>=j ){
            hist_mu_all_count[i]++;
          }
          hist_mu[i]=0;
          hist_count[i]=0;
        }
      /* ANALYSE mu(t) */
      } else {
	if(logt<mu_logtL){
	  printf("Error: value --mu-time-L should not be larger than %e.\n", pow(10, logt)); goto TERMINATE_PROGRAM;
	}
        if(logt>mu_logtU) {
          printf("Error: value --mu-time-U should be at least %e.\n", pow(10, logt)); goto TERMINATE_PROGRAM;
	}
        i=(int)((logt-mu_logtL)/hist_binwidth); /* bin number */
        hist_mu[i]+=mu;
        hist_count[i]++;
      }
    }
  } /* mu(t) histogram updated */
  /*================================================================*/



  /* EXPORT */
  if(ti%Nprint==0) {
    if(opt_export_logt)
      Nprint*=2;
    Rcm=0.0;   /* centre of mass */
    stress=0.0;  /* calculate stress */
    for (i=0; i<n_beads; i++) {
      Rcm+=R[i];
     // if (sequence[i]==1 && state[i]==1)  /* fraction of closed stickers */
     //   tmpf++;
      if(i>0)
        stress+=pow(R[i]-R[i-1], 2)*di; /*stress */
    }
    Rcm/=(n_beads); 
    ///tmpf/=Zs; /* fraction of closed stickers */
    

    ifp=fopen(file_timeprogress, "a"); /* output monomer coordinates relative to the center of mass */
    fprintf(ifp, "%12.5e %12.5e %12.5e %12.5e %12.5e", ttime, R[0]+dcorr-Rcm, R[n_beads-1]-dcorr-Rcm, Rcm, stress);
    if(opt_export_sequence)
      for(i=0; i<n_beads; i++)
        fprintf(ifp, " %12.5e", R[i]-Rcm);
    fprintf(ifp, "\n"); /* end of data line*/
    fclose(ifp);
  }
}  /* ALL TIME STEPS DONE */
/*================================================*/

  printf("p_closed=%f\n", (ttime-ttopen)/ttime);

/* ANALYSE RESULTS */
if(opt_calc_mu) {
  double muII0, muII, dmuI, dmuII;
  for(i=0; i<mu_nbin; i++) {
    if(hist_mu_all[i]>0)
      hist_mu_all[i]/=hist_mu_all_count[i];
  }


  for(i=0; i<mu_nbin; i++) {
    if(hist_mu_all[i]>0) {
      ttime=pow(10, mu_logtL+hist_binwidth*i);
      if (i==0) {
        muII0=(Nx>0?mu_LMcumm[i]/Ncumm:0); /* normalise */
        /*                                    time          muI(t)  dmuI(t)   muII(t)   dmuII(t) */
        printf("%12e %12e %12e %12e %12e\n", ttime, hist_mu_all[i],    0.0,     muII0,      0.0);
      } else{
        dmuI=-4*Ze*pow( ttime, 0.75 )*(hist_mu_all[i] - hist_mu_all[i-1]  )/(   ttime - pow(10, mu_logtL+hist_binwidth*(i-1))  );

        muII=(Nx>0?mu_LMcumm[i]/Ncumm:0);    /* normalise */
        dmuII=-4*Ze*pow( ttime, 0.75 )*(muII - muII0)/(   ttime - pow(10, mu_logtL+hist_binwidth*(i-1))  ); muII0=muII;
        /*                                    time          muI(t)  dmuI(t)   muII(t)   dmuII(t) */
        printf("%12e %12e %12e %12e %12e\n", ttime, hist_mu_all[i],    dmuI,     muII,     dmuII );
      }
    }
  }
   printf("#nrelaxations: %d ; %e\n", count_relaxations, pow(10, mu_logtU));

}
   printf("#open; closed: %d ; %d\n", countOpen, countClose);
  /*===========================================*/
  TERMINATE_PROGRAM: {        /* FREE MEMORY */
    if(                R!=NULL){free(R);}
    if(               dR!=NULL){free(dR);}
    if(             dRdi!=NULL){free(dRdi);}
    if(             dRdt!=NULL){free(dRdt);}
    if(               ks!=NULL){free(ks);}
    if(             state!=NULL){free(state);}
    if(          hist_mu!=NULL){free(hist_mu);}
    if(      hist_mu_all!=NULL){free(hist_mu_all);}
    if(hist_mu_all_count!=NULL){free(hist_mu_all_count);}
    if(       hist_count!=NULL){free(hist_count);}
    if(         sequence!=NULL){free(sequence);}
    if(             outdir!=NULL){free(outdir);}
    if(             file_timeprogress!=NULL){free(file_timeprogress);}
    if(             file_param!=NULL){free(file_param);}
    if(             mu_instant!=NULL){free(mu_instant);}
    if(             mu_LMcumm !=NULL){free(mu_LMcumm);}
    if(             t_alive   !=NULL){free(t_alive);}
  };

  return(0);
}
