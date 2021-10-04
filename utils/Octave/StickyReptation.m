% Evaluate Sticky-Rouse model: G(t), G'(w), G''(w) are known analytically.
% Use two numerical schemes to determine G' and G'' from G(t) and compare
% to analytical result
function StickyReptation()
  close all; clc;
  addpath('StickyRouse') % Tranform G(t) to G'(w) and G''(w)
  addpath('DoubleReptation')  % G(t) for double reptation
  addpath('LikhtmanMcLeish')  % G(t) for reptation - LikhtmanMcLeish
                              % does NOT include Rouse modes
  addpath('getDynamicModuli') % Tranform G(t) to G'(w) and G''(w)


  nFEM=50; % Finite element analysis
  Cnu=0.7;
  tolerance=1e-4;


  Zs=15.0;
  tau_s=0.01;
  Ze=15.0;
  
  tau_d0=3*tau_s*(Zs/2)^2*Ze;     % Sticky-reptation time
  tS    =3*tau_s*(0.5*Zs)*(0.5*Zs); % Relaxation time of sticker strand

  % Time range for G(t)
  trow=10.^(linspace(-log10(tau_s)-5,log10(tau_d0)+1, 50));

  Grep   = Gt_StickyLikhtmanMcLeish(trow/tau_s, Ze, Zs, Cnu, tolerance);
  trow=trow(Grep>0);
  Grep=Grep(Grep>0);
  

  % Sticky Rouse: Relaxation modulus
  Gsr    = Gt_StickyRouse(trow/tS, Ze, Zs );
  Gt  =Gsr+Grep;
  Gsr =Gsr( Gt>0);
  Grep=Grep(Gt>0);
  trow=trow(Gt>0);
  Gt  =Gt(  Gt>0);
  


  % Numerical dynamic moduli via two methods
  %nFEM=length(Gt)/2
  G1G2  =getDynamicModuli_EvansTassieri(trow, Gt, nFEM);
  G1G2_II  =getDynamicModuli_Nobile(trow, Gt, nFEM);



  %===========================================================
  % Plot results
  figure
  subplot(2,1,1)
  %trow=trow(Gt>0); Grep=Grep(Gt>0); Gt  =Gt(Gt>0);
  plot(log10(trow(Gt>0)), log10(Gt(Gt>0)), 'k', 'LineWidth', 2); hold on
 % loglog(trow, Gsr, 'b', 'LineWidth', 2); hold on;
 % loglog(trow, Grep, 'g', 'LineWidth', 2); hold on;
 % legend('G_{tot}', 'G_{sr}') % , 'G_{rep}'
  xlabel('time')
  ylabel('log_{10}G(t)')

  subplot(2,1,2)
  loglog(G1G2(G1G2(:,2)>0,1), G1G2(G1G2(:,2)>0,2), 'r', 'LineWidth', 2); hold on;
  loglog(G1G2(G1G2(:,3)>0,1), G1G2(G1G2(:,3)>0,3), 'g', 'LineWidth', 2); hold on;
  loglog(G1G2_II(G1G2_II(:,2)>0,1), G1G2_II(G1G2_II(:,2)>0,2), 'k', 'LineWidth', 1); hold on;
  loglog(G1G2_II(G1G2_II(:,3)>0,1), G1G2_II(G1G2_II(:,3)>0,3), 'k', 'LineWidth', 1); hold on;
  axis([0.5/tau_d0, 20/tau_s, 1e-2, 2])
  xlabel('frequency')
  ylabel('G`(w), G``(w)')

  %===========================================================

  % Get pre-calculated values g(x) and F(U) 
  tol=1e-4;    % Tolerance
  Nsamples=40; % 

  SET_DESCLOIZEAUX=initialise_desCloizeaux(tol, Nsamples);
  

  % Pre-calculated values
  x=SET_DESCLOIZEAUX.GTABLE(:,1);
  g=SET_DESCLOIZEAUX.GTABLE(:,2);
  U=SET_DESCLOIZEAUX.FTABLE(:,1);
  F=SET_DESCLOIZEAUX.FTABLE(:,2);

  % Interpolated and extrapolated values
  xL=SET_DESCLOIZEAUX.xL;
  xU=SET_DESCLOIZEAUX.xU;
  UL=SET_DESCLOIZEAUX.UL;
  UU=SET_DESCLOIZEAUX.UU;

  xint=10.^linspace(log10(0.1*xL), log10(10*xU), 200);
  Uint=10.^linspace(log10(0.1*UL), log10(10*UU), 200);
  for i=1:200
    gint(i)=g_descloizeaux_interpolate(xint(i), SET_DESCLOIZEAUX);
    Fint(i)=F_descloizeaux_interpolate(Uint(i), SET_DESCLOIZEAUX);
  end


  % Show results
  figure
  subplot(1,2,1)
  loglog(x, g, '.k', 'MarkerSize', 12); hold on
  loglog(xint, gint, 'r'); hold on
  loglog(x, -x+sqrt(x.*(x+sqrt(pi*x)+pi)), 'k');
  legend('pre-calculated', 'interpolated', '-x+sqrt(x.*(x+sqrt(pi*x)+pi))')
  xlabel('x')
  ylabel('g(x)')
  subplot(1,2,2)
  loglog(U, F, '.k', 'MarkerSize', 12); hold on
  loglog(Uint, Fint, 'r'); hold on
  legend('pre-calculated', 'interpolated', 'Location', 'SouthWest')
  xlabel('U')
  ylabel('F(U)')
end


