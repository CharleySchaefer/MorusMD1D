% time in units of the reptation time
% Cost function
% Time input: units of reptation time
function cost=cost_dmu(param, Z, tdata, dmudata,  tolerance)
  Gamma4th=-4.90166680986071058051639321345156210740495699243228244492;
  Gf=param(1); % Renormalisation of elastic modulus
  Tf=param(2); % Renormalisation of reptation time
  Cmu=param(3); % Early-time plateau value (C=1.5+/-0.02 in LM paper)

  C=-3^(1/4)*Cmu/Gamma4th; % (C = 0.4026 = 0.306*3^(1/4) in LM paper)

  Tfi=1./Tf;               % inverse of Tf
  fac_Gf=8*Gf/(pi^2);
  Zi=1./Z;

  % Crossover epsilon value
  p_star=sqrt(Z*0.1);
  dsum=0.0; p=1;
  while p<=p_star
    dsum=dsum+1.0/(p*p);
    p=p+2;
  end
  epsf=(4*C/(1-dsum*fac_Gf))^4;


  dmurow=zeros(size(tdata));
  for i=1:length(tdata)
    ttime=tdata(i);
    
    dsumB=0.0; p=1;
    while p<=p_star
      p2=p*p;
      tmp=exp(-ttime*p2*Tfi);
      dsumB=dsumB+tmp;
      p=p+2;
    end
    dmurow(i)=-dsumB*fac_Gf*Tfi;

    tZ=ttime*Zi; % tZ=time in units of reptation time multiplied by Z,
                 % so the normalisation factor is proportional to Z^4.
    epstZ=epsf*tZ;
    if(epstZ > 15) % capture for small tZ
      gammafac=0;
    else
      gammafac=gamma_fourth(epstZ, tolerance);
    end
    tmp=C*(tZ)^(  0.25);
    dmurow(i) = dmurow(i) + tmp*( ...
                    -epsf*exp(-epstZ)*(epstZ)^(-5/4) ...
                    +gammafac/(4*tZ))*Zi; 
    % Scale the derivative
    dmurow(i)=-4*Z*(tZ)^(0.75)*dmurow(i)/3^0.25;
  end
  cost=dmurow-dmudata;
end


% Get numerical value of Gamma(-1/4, x)
function Gamma = gamma_fourth(x, tolerance)
  Gamma4th=-4.90166680986071058051639321345156210740495699243228244492;
  fac=x^(-0.25);
  kfac=1.0;
  dsum=-4.0*fac;
  product=1.0*fac;
  B=product/kfac;
  Gamma=Gamma4th-dsum;
  err=2*tolerance;
  k=0;
  while (err>tolerance)
    k=k+1;
    %kfac=kfac*k;
    %product=product*-x;
    B=B*(-x/k); %product/kfac

    ddsum=B/( (k-0.25));
    Gamma=Gamma-ddsum;
    err=abs(ddsum/Gamma);
  end
end
