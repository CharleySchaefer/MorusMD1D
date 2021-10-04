% Evaluate Sticky-Rouse model: G(t), G'(w), G''(w) are known analytically.
% Use two numerical schemes to determine G' and G'' from G(t) and compare
% to analytical result
function demo_getModuli()
  close all; clc;
  addpath('StickyRouse') % Tranform G(t) to G'(w) and G''(w)
  addpath('LikhtmanMcLeish')  % G(t) for reptation - LikhtmanMcLeish
                              % does NOT include Rouse modes
  addpath('DoubleReptation')  % 
  addpath('getDynamicModuli') % Tranform G(t) to G'(w) and G''(w)

  %==============================================================
  % USER SETTINGS
  % A. Polymer properties
  Ze=100;       % Number of entanglements per chain
  tau_e=1e-6;  % Rouse time of entanglement strand
               % NOTE: tau_e not used if Zs>0 
  Zs=0;        % Number of stickers per chain
  tau_s=1e-2;  % Sticker dissociation time

  % B. Materials parameters
      % Likhtman-McLeish constitutive model
  LM_Cnu=0.010;  % Likhtman-McLeish: strength of contour-length fluctuations
      % desCloizeaux     - Double-Reptation model
  DR_alpha=4.0;
  DR_beta=2.00;

  % C. Experimental settings
  tL = 1e-4; % Shortest time (1/highest frequency)
  tU = 1e4;  % Longest time  1/lowest frequency) 
  Nt = 500;   % Number of data points

  % D. Numerical parameter values
  nFEM=40;           % Number of frequencies to calculate G(w)
  LM_tolerance=1e-5; % Numerical tolerance for LikhtmanMcLeish
  initialise_DR=0; % 0: load pre-calculated  Double Reptation
                   % 1: new  pre-calculation Double Reptation; 
  % DR precalculation parameters 
  Nsamples=40; % Number of samples: g(x) and F(U)
  tol=1e-4;    % Numerical tolerance for g(x) and F(U) calculation
  Nt=64;       % Number of samples time 
  NH=32;       % Number of samples H parameter
  DR_HL=1e-2;  DR_HU=1e2; % Interpolation interval: H
  DR_tL=1e-6;  DR_tU=1e2; % Interpolation interval: t
  % END DR precalculation parameters 
  %==============================================================

  if initialise_DR==0 % Load precalculated DR model
    load('DoubleReptation/Precalculations');
  elseif initialise_DR==1
    DOUBLE_REPTATION = initialise_DoubleReptation( Nsamples, Nt, NH, tol, DR_HL, DR_HU, DR_tL, DR_tU);
    save('DoubleReptation/Precalculations', 'DOUBLE_REPTATION');
  end

  % Time range for G(t)
  trow=10.^(linspace(log10(tL),log10(tU), Nt));

  % Physical time scales
  if Zs==0
    tau_d0=3*tau_e*Ze^3;    % Reptation time
  elseif (Zs>0)
    tau_d0=3*tau_s*Zs^2*Ze; % Sticky-reptation time
    tS=tau_s*(Zs)*(Zs);     % Relaxation time of sticker strand
  else
    error('Zs cannot be negative.\n');
  end
  fprintf('reptation time: %e\n', tau_d0);


  % REPTATION
  Grep=zeros(size(trow));  % In case reptation is switched off

  % Double reptation
  H    = (Ze/DR_alpha);
  tt   = (trow/tau_d0);
  GDR  =10.^((DR_beta/2)*Grep_desCloizeaux_interpolate(tt, H, DOUBLE_REPTATION));
  G1G2DR  =getDynamicModuli_EvansTassieri(trow, GDR, nFEM);

  % Likhtman-McLeish reptation
  GLM = Gt_LikhtmanMcLeish(trow/tau_e, Ze, LM_Cnu, LM_tolerance);
  G1G2LM  =getDynamicModuli_EvansTassieri(trow, GLM, nFEM);
  
  % derivative
  trow_d=zeros(1,length(trow)-1);
  GDR_d =zeros(1,length(trow)-1);
  GLM_d =zeros(1,length(trow)-1);
  for i =1:length(trow)-1
    trow_d(i)=trow(i+1);
    dt=trow(i+1)-trow(i);
    GDR_d(i)=-4*Ze*tau_e^0.25*trow(i+1)^0.75*( GDR(i+1)-GDR(i) )/dt;
    GLM_d(i)=-4*Ze*tau_e^0.25*trow(i+1)^0.75*( GLM(i+1)-GLM(i) )/dt;
  end


  %===========================================================
  % Plot results
  figure
  subplot(2,2,1)
  loglog((trow(GDR>0)), (GDR(GDR>0)), 'k', 'LineWidth', 2); hold on
  loglog((trow(GLM>0)), (GLM(GLM>0)), '--k', 'LineWidth', 2); hold on
  axis([tau_e, 2*tau_d0, 1e-5, 2])
  legend(sprintf('Double Reptation; alpha=%.1f, beta=%.2f', DR_alpha, DR_beta), sprintf('Likhtman-McLeish; Cnu=%.2f', LM_Cnu), 'Location', 'SouthWest')
  xlabel('time')
  ylabel('log_{10}G(t)')

  subplot(2,2,2)
  plot( log10(trow_d(GDR_d>0)), GDR_d(GDR_d>0), 'k', 'LineWidth', 2 ); hold on;
  plot( log10(trow_d(GLM_d>0)), GLM_d(GLM_d>0), '--k', 'LineWidth', 2 )
  xlabel('time')
  ylabel('dG')

  %loglog(trow, Gsr, 'b', 'LineWidth', 2); hold on;
  %loglog(trow, Grep, 'g', 'LineWidth', 2); hold on;
 % legend('G_{tot}', 'G_{sr}') % , 'G_{rep}'


  subplot(2,2,3:4)
  loglog(G1G2DR(G1G2DR(:,2)>0,1), G1G2DR(G1G2DR(:,2)>0,2), 'r', 'LineWidth', 2); hold on;
  loglog(G1G2DR(G1G2DR(:,3)>0,1), G1G2DR(G1G2DR(:,3)>0,3), 'g', 'LineWidth', 2); hold on;
  loglog(G1G2LM(G1G2LM(:,2)>0,1), G1G2LM(G1G2LM(:,2)>0,2), '--r', 'LineWidth', 2); hold on;
  loglog(G1G2LM(G1G2LM(:,3)>0,1), G1G2LM(G1G2LM(:,3)>0,3), '--g', 'LineWidth', 2); hold on;
  %axis([0.5/tau_d0, 20/tau_s, 1e-7, 2])
  xlabel('frequency')
  ylabel('G`(w), G``(w)')

 
end

