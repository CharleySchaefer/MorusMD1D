%Rubinstein&Colby - JCP 89, 5291 (1988), https://doi.org/10.1063/1.455620

% The contrain release rate can be calculated as R(t) = \int deps * dM(eps)/deps *exp(-eps * t).
% This script calculated M(eps) from mu(t): First a distribution P(eps)deps is obtained through a Laplace transform and M(eps) is calculated from P(eps)d(eps) through the algorithm in Appendix C, with M(eps) in Fig 8 as a result.
function SCT_Rubinstein()
  close all; clc;

  %----------------
  % USER INPUT
  N=200;
  nu=1.0
  tau=1.0;
  %----------------

  % PRECALCULATIONS
  eps_star1=1.0/(tau*(1-nu/sqrt(N))^2)
  eps_star2=N/nu^2/tau
  prefac1=0.5/sqrt(tau);  
  prefac2=0.5*(nu^2/(N*tau))^0.25;  
  sum_cont1=prefac1*(2/sqrt(eps_star1)-2/sqrt(eps_star2) ) % Integral
  sum_cont2=4*prefac2*eps_star2^(-0.25)  % Integral


  % CUMMULATIVE PROBABILITY
  eps_row=10.^linspace(log10(min(eps_star1))-4, log10(max(eps_star2))+9, 100);
  S=zeros(1,100);
  for i=1:100
    if eps_row(i)<=eps_star1
      P(i) = 0;
      S(i)=S(i)+P(i);
      P2(i) = 0;
    elseif eps_row(i)<=eps_star2
      S(i)=2*prefac1*(eps_star1^(-0.5)-eps_row(i)^(-0.5));
      P(i) = prefac1*eps_row(i)^(-1.5)*(eps_row(i)-eps_row(i-1));
    else
      S(i)=sum_cont1 + 4*prefac2*(eps_star2^(-0.25)-eps_row(i)^(-0.25));
      P(i)  = prefac2*eps_row(i)^(-5.0/4)*(eps_row(i)-eps_row(i-1));
    end
  end

  % SAMPLING FUNCTION
  P_row=[0.01:0.01:1];
  e_row=zeros(size(P_row));
  for i=1:length(P_row)
    if P_row(i)<sum_cont1
      e_row(i) = ( eps_star1^(-0.5) - P_row(i)/(2*prefac1) )^-2;
    else 
      e_row(i) = ( eps_star2^(-0.25) - (P_row(i) - sum_cont1)/(4*prefac2))^-4;
    end
  end 


  % SAMPLE
  Nsamples=4000;
  eps_u=zeros(1,Nsamples);sum1=0;sum2=0;
  rnd=rand(1,Nsamples);
  test2 =   sum(rnd<sum_cont1)/Nsamples
  for j=1:Nsamples
    u=rand()  ;
    if u<sum_cont1
      eps_u(j) = ( eps_star1^(-0.5) - u/(2*prefac1) )^-2;
    else 
      eps_u(j) = ( eps_star2^(-0.25) - (u - sum_cont1)/(4*prefac2))^-4;
    end
  end

  % GET M(eps) FUNCTION
  s=zeros(1,N-1);
  M=zeros(length(eps_row));
  N_chains=1;
  for r=1:N_chains
    m=sample_chain_mobilities(N, sum_cont1, eps_star1, prefac1, eps_star2, prefac2);
  for k=1:length(eps_row)
    s(1)=m(1)+m(2)-eps_row(k);
    for i=2:N-1
      s(i)=m(i)+m(i+1)-eps_row(k)-m(i)^2/s(i-1);
    end
    M(k)=M(k)+sum(s<0);
  end
  end
  M=M/N_chains;

  figure
  subplot(2,2,1)
    plot(log(eps_row(S>0)), S(S>0))
    xlabel('\epsilon')
    ylabel('S')
  subplot(2,2,2)
    plot(P_row(P_row>0), log(e_row(P_row>0)) )
    ylabel('\epsilon')
    xlabel('S')
  subplot(2,2,3)
    hist(log10(eps_u), 200); hold on
 %   plot(log(eps_row(count>0)), count(count>0), '.k'); hold on
    plot(log10(eps_row(P>0)), P(P>0)*Nsamples/1.5, 'r'); hold on
  xlabel('\epsilon')
  ylabel('P(\epsilon)d(\epsilon)')
  set(gca,'YScale','log')
  subplot(2,2,4)
  plot(log10(eps_row*sqrt(N)), M/N)


end

function m=sample_chain_mobilities(N, sum_cont1, eps_star1, prefac1, eps_star2, prefac2)
  u=rand(1,N); m =zeros(size(u));
  for j=1:N
    if u<sum_cont1
      m(j) = ( eps_star1^(-0.5) - u(j)/(2*prefac1) ).^-2;
    else 
      m(j) = ( eps_star2^(-0.25) - (u(j) - sum_cont1)/(4*prefac2)).^-4;
    end
  end
end

