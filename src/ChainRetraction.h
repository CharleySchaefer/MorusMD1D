/* Pade approximant for the inverse Langevin function [Cohen, Rheol Acta (1991)] */
double get_ks_norm(double lambda_max) {
  double lmax2=lambda_max*lambda_max;
  return (lmax2-1.0)/(3*lmax2-1);   /* ks_norm */
}

/*

  Analytical expression is given by

    ks(lambda) = norm_fac*( 3*lambda_max^2 - lambda^2 )/( lambda_max^2 - lambda^2 )
  
    here, norm_fac is defined such that 

    ks(lambda=1) = 1

  To numerically deal with the singularity at lambda = lambda_max,
  we only use the analytical expression up to a cutoff at 'l_cutoff'< lambda_max;

  for lambda>l_cutoff, we use a an exponential extrapolation

    ks(lambda) = A*exp( B*(lambda-lambda_cutoff) ),

  with 

    A = norm_fac*( 3*lambda_max^2 - lambda_cutoff^2 )/( lambda_max^2 - lambda_cutoff^2 )

  for continuity of the function ks(lambda),   and

    B = norm_fac*(4 lambda_cut lambda_max^2)/(     lambda_cut^4 -  4 lambda_cut^2 lambda_max^2 + 3 lambda_max^4     )

  for continuity of the derivative ks'(lambda)


*/
double get_ks(double lambda, double lambda_max, double ks_norm){
  double ks, A,B;
  double l_cutoff=0.99*lambda_max, l_cutoff2=l_cutoff*l_cutoff;
  double lmax2=lambda_max*lambda_max, denomi;
  double l2=lambda*lambda;

  if(lambda>l_cutoff) { 
    A=ks_norm*( 3*lmax2 - l_cutoff2 )/( lmax2 - l_cutoff2 );
    B=ks_norm*( 4*l_cutoff*lmax2 )/(     l_cutoff2*(l_cutoff2 -  4*lmax2) + 3*lmax2*lmax2     );
    ks=A*exp( B*(lambda-l_cutoff) );

  } else {
    ks=ks_norm*(3.0*lmax2-l2)/(lmax2-l2);  /* ks_norm */
  }
  return ks;
}
