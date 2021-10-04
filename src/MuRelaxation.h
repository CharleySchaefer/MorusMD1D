/*
  Reference: Likhtman-McLeish paper
*/

/* 
  discretise space to track lifetime of alive chain segments
    R:     x-coordinate in (-inf, inf)
    x[i]:  periodic in [0 to Lx),
           with grid spacing dx.
           x[i] = dx*i
*/
int get_grid_index( double R, int Nx, double dx ){
  int xi=((int)round(  R/dx ))%Nx; /* discretise */
  if( xi<0 )                /* shift to positive value */
    xi+=Nx;
  return xi;
}



/*
  mu_instant( t ), with t discretised time on log scale:

    t[i]=10^( logtL + i*dlogt )

  For a given time, the index i is found as

      i = (int)( ( log10(time) - logtL )/dlogt )
*/
int get_mu_instant(double *mu_instant, int Nt, double logtL, double dlogt, double *t_alive, int Nx, int xiL, int xiU ) {

  int i,j,k, 
      imax=-1; /*output -1 in case of error */
  int Nalive, loop;
  double logt;
  /*initialise */
  for(i=0; i<Nt; i++) {
    mu_instant[i]=0.0;
  }

  Nalive=0;                           /* count number of alive segments       */
  j=xiL;                              /* first segment */
  /* loop over alive segments on a periodic grid */
  while(1) {                          /* will be ended using "break;"         */
    /* analyse grid segment j */
    Nalive++;                         
    logt=log10(t_alive[j]  );         /* lifetime of grid segment on log scale*/
    if( logt>=logtL ){
      i = (int)( (logt-logtL)/dlogt  ); /* time index                         */
      imax=(i>imax? i:imax);            /* keep track of longest lifetime     */
      for ( k=0; k<=i; k++ )            /* count alive segment                */
        mu_instant[k]++;
    }
    /* next grid segment */
    if(j==xiU) {                          /* last segment analysed - done */
      break;    /* exit while loop */
    } else {
      j++;
      if(j==Nx) /* periodic boundary */
        j=0;
    }
  } /*end while loop */ 

  for ( k=0; k<=imax; k++ )           /* normalise mu(t)                      */
    mu_instant[k]/=(Nalive); 

  return imax;
}
