% Calculate Rate of Constraint Release, R(t), for Likhtman-McLeish theory.
% This script reproduces Fig 6 in the Likhtman-McLeish paper.
% The algorithm is in more detail discussed in the Rubinstein-Colby 1988 paper, see SCT_Rubinstein.m

function SCT_LikhtmanMcLeish()
  close all; clc;
  Cnu =1.0;
  Z   =300;
  taue=1.0;

  taud0=3*taue*Z^3; % Reptation time (without CLF)

  % CLF
  Gf=get_renormalisation_elastic_modulus(Z);  
  Tf=get_renormalisation_reptation_time(Z);   
  taudf=taud0*Tf; 
  
  % Crossover from discrete to continuous disctribution
  epsf=get_eps_star(Z,Gf); % in units of taud0*Z
  epsf=epsf/(taud0*Z)      % Cutoff value

  % Discrete contribution (Doi-Edwards delta-function peaks)
  prefac1=8*Gf/pi^2;
  m=1; i=1; cummP=zeros(1,sqrt(Z/10));
  while m <= sqrt(Z/10)
    m2=m*m;
    epsilon(i) = m2/taudf;
    P(i)       = prefac1/m2;
    if i==1
      cummP(i)=P(i);
    else
      cummP(i)=cummP(i-1)+P(i);
    end
    i=i+1; 
    m=m+2; % odd m
  end
  cummP(i:end)=[];
  sum_discreteP=sum(P); % Integral over discrete contributions

  % integral continuous part:
  prefac2=0.4026/((Z*taud0)^0.25)  
  sum_contP=4*prefac2*epsf^(-0.25)  % Integral over continuous part
  sumP=sum_contP+sum_discreteP;     % Check consistency

  % Discretise continuous part
  epsilon_continuous=10.^linspace(log10(epsf), log10(epsf)+5, 100);
  eps_c_d=0.5*(epsilon_continuous(1:end-1)+epsilon_continuous(2:end));
  d_eps_c=-(epsilon_continuous(1:end-1)-epsilon_continuous(2:end));
  P_continuous = prefac2*epsilon_continuous(2:end).^(-5/4).*d_eps_c;
  cummP_c = sum_discreteP + 4*prefac2*( epsf^(-0.25)- eps_c_d.^(-0.25) );

  %===============================
  % STEP 1 GET P(E)
  % draw uniform random number between 0 and 1
  fprintf('TEST sampling of P(eps).\n')
  tic
  Nsamples=4000;
  for j=1:Nsamples
    u=rand()  ;
    if u <= sum_discreteP % Dirac-delta
      i=1;
      while ( cummP(i)<u) % i<length(cummP) &&
        i=i+1;
      end
        eps_u(j)=epsilon(i);
    else % continuous part
      
 % cummP_c = sum_discreteP + 4*prefac2*( epsf^(-0.25)- eps_c_d.^(-0.25) );
        eps_u(j)= (epsf^(-0.25) - (u-sum_discreteP)/(4*prefac2))^-4;
    end
  end
  P_cumm_cont=4*prefac2*epsilon_continuous.^-0.25;
  toc

  figure
  subplot(2,2,1:2)
  loglog(epsilon(1:length(P)),            P,           '.k','MarkerSize', 12); hold on;
  loglog(eps_c_d, P_continuous,'.r','LineWidth', 2)
  xlabel('\epsilon')
  ylabel('P(\epsilon)d\epsilon')
  
  subplot(2,2,3)


  plot(log10(epsilon(1:length(cummP))),  cummP, '.k'); hold on;
 % subplot(2,2,4)
  plot(log10(eps_c_d), cummP_c, 'r');
  %loglog(P_cumm_cont,  (P_cumm_cont/(4*prefac2)).^-4, 'k')

  figure
  Nbins=100;
  hist(log10(eps_u), Nbins); hold on;
  plot(log10(epsilon_continuous(P_continuous>0)), (P_continuous(P_continuous>0)), 'k')
 % set(gca,'YScale','log')
  xlabel('epsilon')
  ylabel('P(epsilon)')

  %=============================
  % GENERATE CHAIN

  fprintf('Solve (stochastic) recurvsive equation.\n');
  tic

  Nsamples=300;
  Neps=50;     
  eps_row=10.^linspace(-10, 0, Neps);   

  mrow=zeros(1,Z+1); % Z random numbers


  for j=1:Z+1
    u=rand()  ;
    if u <= sum_discreteP % Dirac-delta
      i=1;
      while ( cummP(i)<u) % i<length(cummP) &&
        i=i+1;
      end
        mrow(j)=epsilon(i);
    else % continuous part
        mrow(j)= (epsf^(-0.25) - (u-sum_discreteP)/(4*prefac2))^-4;
    end
  end





  % STEP 2 - Get M(epsilon)
  % Computational time scales linear with Nsamples and Neps.

 
    mrow=zeros(1,Z+1); % Z random numbers
  Mcumm=zeros(1,Neps);
  %warning('off'); % Switch of warnings (occur due to division by zero)
  for l=1:Nsamples
    % Get random chain mobility
    for j=1:Z+1
      u=rand()  ;
      if u <= sum_discreteP % Dirac-delta
        i=1;
        while ( cummP(i)<u) % i<length(cummP) &&
          i=i+1;
        end
          mrow(j)=epsilon(i);
      else % continuous part
          mrow(j)= (epsf^(-0.25) - (u-sum_discreteP)/(4*prefac2))^-4;
      end
    end


    % Loop eps values ( much slower than get_chain_mobility() )
    sinf=0;
    for k=1:Neps % Loop eps values

      % Initial values recursive relation:  s(1)= m(1)+m(2)-e
      mj=mrow(1);                            % m(1)
      mjp=mrow(2  );                         % m(2)
      epsk=eps_row(k);                       % e
      sprev=mj+mjp-epsk;                     % s(1) 
      if(sprev<0)    % collect statistics (count negative s(j)'s)
        count=1;
      else
        count=0;
      end;

      % Loop recursive relation        % s(j)=m(j)+m(j+1)-e -mj^2/s(j-1)
      for j=2:1:Z
        mj=mjp;                          % m(j)
        mjp=mrow(j+1);                   % m(j+1)

        if (sprev==0) % Capture division by 0 
          sinf=1; % infinity encountered
          count=count+1; % result -inf --> negative s(j) is counted
        elseif (sinf==0)
          sprev=mj+ mjp - epsk - mj*mj/sprev;
          if(sprev<0)  % collect statistics (count negative s(j)'s)
            count=count+1;
          end;
        else
          sinf=0; % reset
          sprev=mj+ mjp - epsk; 
          if(sprev<0)  % collect statistics (count negative s(j)'s)
            count=count+1;
          end;
        end
      end

      % Save statistics
      Mcumm(k)=Mcumm(k)+count;
    end
  end
  %warning('on') % Switch warnings back on
  toc

  %M=M/Z;
  %M=mean(M) %Mstd=zeros(size(Mmean));
  M=Mcumm/(Nsamples*Z);
  %for l=1:Nsamples
  %  Mstd=(M(l,:)-Mmean).^2;
  %end
  %Mstd=sqrt(mean(Mstd))

  % STEP 3 - Get R(t)
  dMde=zeros(1,Neps-1);
  deps=zeros(1,Neps-1);
  eps_d=zeros(1,Neps-1);
  %eps_d=(eps_row(1:Neps-1)+eps_row(2:Neps))*0.5;
  for i=1:Neps-1
    eps_d(i)= ( eps_row(i+1) + eps_row(i) )/2;
    deps(i) = ( eps_row(i+1) - eps_row(i) );
    dMde(i) = (       M(i+1) -       M(i) )/deps(i); %/deps(i);%; %*;
  end

  % Integrate
  Ntime=50;
  trow=10.^linspace(0, 10, Ntime); % units of tauE
  for i=1:Ntime
    R(i)=  sum( dMde.*exp(-Cnu*eps_d*trow(i)).*deps ); %;
  end

  % mu(t)
  tolerance=1e-4;
  murow=mu_LikhtmanMcLeish(trow/taud0, Z, tolerance); 

  % Differentiate
  dRdt =zeros(1,Ntime-1);
  dmudt =zeros(1,Ntime-1);
  dtrow=zeros(1,Ntime-1); 
  trow_d=zeros(1,Ntime-1); 
  for i=1:Ntime-1
    trow_d(i) = (trow(i+1)+trow(i))/2;
    dtrow( i) = (trow(i+1)-trow(i));
    dRdt(  i) = (R(i+1)-R(i))/dtrow(i);
    dmudt(  i) = (murow(i+1)-murow(i))/dtrow(i);
  end

  Rearly=1-1.8/Z*(Cnu*trow).^0.25; % LM 

  dRdt_scaled =-4*Z*taue^0.25*trow_d.^0.75.*dRdt;
  dmudt_scaled=-4*Z*taue^0.25*trow_d.^0.75.*dmudt;
  figure
  subplot(2,2,1)
  plot(log10(eps_row),M); hold on
 % plot(log10(eps_row),1-exp(-(eps_row*10000).^0.5));
  axis([-10,4,0,1])
  xlabel('\epsilon')
  ylabel('M')
  subplot(2,2,2)
  plot(log10(eps_d),dMde);
  xlabel('epsilon')
  ylabel('dM/d\epsilon')
  subplot(2,2,3)
  loglog(trow(R>1e-8),R(R>1e-8)); hold on
  loglog(trow(Rearly>1e-8),Rearly(Rearly>1e-8), '--k');
  legend('numerical', 'approx', 'Location', 'southwest')
  xlabel('time')
  ylabel('R(t)')
  subplot(2,2,4)
  loglog(trow_d( dRdt_scaled>1e-8),dRdt_scaled(  dRdt_scaled>1e-8));hold on;
  loglog(trow_d(dmudt_scaled>1e-8),dmudt_scaled(dmudt_scaled>1e-8));
  xlabel('time')
  ylabel('dRdt')
  axis([1e0,1e10,0.1,4])
  
end


