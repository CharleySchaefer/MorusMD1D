% Calculate Relaxation by Constraint Release, R(t), 
% for Likhtman-McLeish theory.
% Algorithm tested in SCT_*
% input time: units of reptation time
function R = R_LikhtmanMcLeish(trow, Z, Nbeads)
  %--------------------------------------------------------
  % USER INPUT
  % Physical Parameters
  % Cnu =1.0; % Cnu just shifts the solution on the x-axis
  %Z      = 30; % Number of entanglements per chain
  %trow;        % Time range of result (units of taud0)

  % Numerical settings
  %Nbeads = 30; % Number of beads to discretise Z  
  %     Discretise epsilon space
  Neps=50;     
  epsMin=-2; % shifted w.r.t. eps_star
  epsMax=2;
  %     Solve stochastic equation for M(epsilon)
  NsamplesM=300;
  %---------------------------------------------------------


  %===============================
  % STEP 1 GET P(epsilon)
  taud0=3*Z^3; % Reptation time (without CLF)
  Gf=get_renormalisation_elastic_modulus(Z);  
  Tf=get_renormalisation_reptation_time(Z);   
  taudf=taud0*Tf; 
  
  % Crossover from discrete to continuous disctribution
  epsf=get_eps_star(Z,Gf); % in units of taud0*Z
  epsf=epsf/(taud0*Z);      % Cutoff value

  % Discrete contribution (Doi-Edwards delta-function peaks)
  prefac1=8*Gf/pi^2;
  m=1; i=1; cummP=zeros(1,sqrt(Z/10));
  while m <= sqrt(0.1*Z)
    m2=m*m;
    eps_discrete(i) = m2/taudf;
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
  Pdiscrete=cummP(end);

  eps_row=10.^linspace(epsMin+log10(min(eps_discrete)), epsMax, Neps);  

  %=============================
  % STEP 2 - Get M(epsilon)
  % Computational time scales linear with Nsamples and Neps.
  prefac2=0.4026/((Z*taud0)^0.25)  ;
  mrow=zeros(1,Nbeads+1); % Z random numbers
  M=zeros(1,Neps);
  for l=1:NsamplesM
    % Get random chain mobility
    for j=1:Nbeads+1
      u=rand()  ;
      if u <= Pdiscrete % Dirac-delta
        i=1;
        while ( cummP(i)<u) % i<length(cummP) &&
          i=i+1;
        end
          mrow(j)=eps_discrete(i);
      else % continuous part
          mrow(j)= (epsf^(-0.25) - (u-Pdiscrete)/(4*prefac2))^-4;
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
      for j=2:1:Nbeads
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
      M(k)=M(k)+count;
    end
  end
  M=M/(NsamplesM*Nbeads);

  %=============================
  % STEP 3 - Get R(t)
  dMde=zeros(1,Neps-1);
  deps=zeros(1,Neps-1);
  eps_d=zeros(1,Neps-1);
  for i=1:Neps-1
    eps_d(i)= ( eps_row(i+1) + eps_row(i) )*0.5;
    deps(i) = ( eps_row(i+1) - eps_row(i) );
    dMde(i) = (       M(i+1) -       M(i) )/deps(i);
  end

  % Integrate -> Final result
  Ntime=length(trow)
  R=zeros(1,Ntime);
  for i=1:Ntime
    R(i)=   (dMde.*exp(  -eps_d*trow(i)*taud0  ))*deps' ; 
    % integral over epsilon
  end 
end


