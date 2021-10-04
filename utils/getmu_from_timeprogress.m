function getmu()
  close all; clc;
  tolerance=1e-7;
  Z=20.0;
  log_export=0;
  Nsteps=30000;
  data=importdata('ttestD.txt');
  data=importdata('test-fixed.out');
  t=data(:,1);
  RL=data(:,2);
  RU=data(:,3);
  dt=t(2)-t(1);
  Ndata=length(t);
  [t_all, mu_all]=mu_simulation(RL, RU, Z, Nsteps, log_export);

  Nsteps3=400;
  Z3=100.0;
  data=importdata('ttest.txt');
  t=data(:,1);
  RL=data(:,2);
  RU=data(:,3);
  Ndata=length(t);
  dt3=t(2)-t(1);
  [t_all3, mu_all3]=mu_simulation(RL, RU, Z3, Nsteps3, log_export);

  tth=10.^[-5:0.05:7.75];
  [muth]=mu_theory(tth, 1.0, Z, tolerance);
  tth3=10.^[-5:0.05:7.75];
  [muth3]=mu_theory(tth3, 1.0, Z3, tolerance);

  figure
  plot(log10(t_all*dt), mu_all, '.g', 'MarkerSize', 10); hold on
  plot(log10(t_all3*dt3), mu_all3, '.k', 'MarkerSize', 10); hold on
  plot(log10(tth), muth, 'r', 'LineWidth', 2); hold on
  plot(log10(tth3), muth3, 'r', 'LineWidth', 2); hold on
  %axis([min(log10(dt*delta_row)),max(log10(dt*delta_row)),0,1])
end

function mu = mu12(RL, RU, t1, t2)
  if(      RL(t1)< RL(t2) && RL(t2)<RU(t1) )
    mu = RU(t1)-RL(t2);
  elseif ( RL(t1)< RU(t2) && RU(t2)<RU(t1) )
    mu = RU(t2)-RL(t1);
  else
    mu=0.0;
  end
end
function [t_all, mu_all]=mu_simulation(RL, RU, Z, Nsteps, log_export)
  Ndata=length(RL)
  t1=1;
  delta=1;
  t_all=[t1]; t0=0;
  while t1-t0<Nsteps 
      if log_export==0
        delta=round(1+delta*1.5);
      end
      t1=t0+delta;
    t_all=[t_all, t1];
  end
  mu_all=zeros(size(t_all));

  Niter=500;
  for j=1:Niter;

    t0=randi([1,Ndata-Nsteps],1,1);
    maxRL=RL(t0);
    minRU=RU(t0);
    delta=1;
    mu_all(1)=1;t1=t0+1; i=1;
    delta=1; mu=1;
    while mu>0 && t1-t0<Nsteps 
      maxRL=max(maxRL, RL(t1));
      minRU=min(minRU, RU(t1));
      mu=(minRU-maxRL);
      i=i+1;
      mu_all(i)=mu_all(i)+mu;
      if log_export==0
        delta=round(1+delta*1.5);
      end
      t1=t0+delta;
    end
  end
  mu_all=mu_all/(Z*Niter);
end

function [ murow]=mu_theory(trow, tau_e, Z, tolerance)

  % Solve Eq.(14)
  p_star=sqrt(Z*0.1);
  dsum=0.0; p=1;
  while p<p_star
    dsum=dsum+1.0/(p*p);
    p=p+2;
  end
  Gf=get_Gf(Z);
  FAC=8*Gf/(pi*pi);
  eps_star=(4*0.306/(Z*(1-FAC*dsum)))^4/tau_e
  
  
  % Solve Eq.(13)
  taudf=get_taudf(Z, tau_e);

  %ttime=1.0e1
  
  for i=1:length(trow)
    ttime=trow(i);
    
      dsum=0.0; p=1;
  while p<p_star
    dsum=dsum+exp(-ttime*p*p/taudf)/(p*p);
    p=p+2;
  end
  
  mu=0*FAC*dsum +0.306/(Z*tau_e^0.25)*ttime^0.25*gamma_fourth(eps_star*ttime, tolerance);
    murow(i)=mu;
  end
end

function taudf=get_taudf(Z, tau_e)
    sqZi=1.0/sqrt(Z);
  taudf=3*tau_e*Z*Z*Z*(1-2*1.69*sqZi+4.17*sqZi*sqZi-1.55*sqZi*sqZi*sqZi);
end

function Gf=get_Gf(Z)
  sqZi=1.0/sqrt(Z);
  Gf=1.0+(-1.69 + (2.0 - 1.24*sqZi)*sqZi)*sqZi;
end

function Gamma = gamma_fourth(tau, tolerance)
  if(tau>12)
    Gamma=0;
  else
  Gamma4th=-4.90166680986071058051639321345156210740495699243228244492;
  fac=tau^(-0.25);
  kfac=1;
  dsum=-4.0*fac;
  product=1.0*fac;
  Gamma=Gamma4th-dsum;
  err=2*tolerance;
  k=0;
  while (err>tolerance)
    k=k+1;
    kfac=kfac*k;
    product=product*-tau;
    ddsum=product/( (k-0.25)*kfac);
    Gamma=Gamma-ddsum;
    err=abs(ddsum/Gamma);
  end
  end
end
