function TwoStateStickyDEMG()

  %----------------
  % Chain properties
  tauS=10000;
  ZE=10;
  p=0.9;
  % Flow rate
  doteps=0.0001; % without stickers (p=0), doteps=0.01 --> Wi=doteps*ZE*ZE=1
  %---------------

  tL=1; % units of tauE
  tU=10*tauS;
  Nt=100;
  trow=10.^linspace( log10(tL), log10(tU), Nt );

  % Solve Master Equation
  P1=getP(trow, 0, tauS, p);
  P2=getP(trow, 1, tauS, p);

  % Solve DEMG
  tauR=ZE*ZE; % Rouse time
  Nt=1000;
  dt=10.0;
  tt=zeros(1,Nt);
  l1=ones(1,Nt);
  l2=ones(1,Nt);

  for i=2:Nt
    tt(i)=tt(i-1)+dt;
    dl1dt=doteps*l1(i-1) + getP(tt(i), 0, tauS, p)/tauR*(1-l1(i-1));
    dl2dt=doteps*l2(i-1) + getP(tt(i), 1, tauS, p)/tauR*(1-l2(i-1));
    l1(i)=l1(i-1)+ dt*dl1dt;
    l2(i)=l2(i-1)+ dt*dl2dt;
  end

  ifp=fopen(sprintf('p%fe%f.out', p, doteps), 'w');
  for i=1:Nt
    fprintf(ifp, '%12e %12e %12e %12e\n', tt(i), l1(i), l2(i), p*l1(i)+(1-p)*l2(i));
  end
  fclose(ifp);
end

% P:    Probability that a sticker is open at time t
% t:    time in units of tauE
% tauS: sticker life time in units of tauE
% p:    equilibrium fraction of closed stickers
function P=getP( t, P0, tauS, p )
  kopen  =1.0/tauS;     % Rate of sticker opening
  kclose=kopen*p/(1-p); % Rate of sticker closing
  k=kopen+kclose;
  K=kopen/k;

  P = K - ( K-P0  )*exp(-k*t);
end
