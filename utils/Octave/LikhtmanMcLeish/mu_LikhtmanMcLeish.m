% time in units of the reptation time
function [murow, dmurow]=mu_LikhtmanMcLeish(trow, Z, tolerance)

  % Influence of contour-length fluctuations
  Gf  =get_renormalisation_elastic_modulus( Z); 
  Tf  =get_renormalisation_reptation_time(  Z); 
  Tfi=1./Tf;               % inverse of Tf
  epsf=get_eps_star(Z,Gf); % units of reptation time*Z
  p_star=sqrt(0.1*Z);
  fac_Gf=8*Gf/(pi^2);
  Zi=1./Z;
  
  for i=1:length(trow)
    ttime=trow(i);
    
    dsum=0.0;dsumB=0.0; p=1;
    while p<=p_star
      p2=p*p;
      tmp=exp(-ttime*p2*Tfi);
      dsum =dsum +tmp/p2;
      dsumB=dsumB+tmp;
      p=p+2;
    end
    
    murow( i)=  dsum*fac_Gf;
    dmurow(i)=-dsumB*fac_Gf*Tfi;

    tZ=ttime*Zi; % tZ=time in units of reptation time multiplied by Z,
                %    so the normalisation factor is proportional to Z^4.
    epstZ=epsf*tZ;
    if(epstZ > 15) % capture for small tZ
      gammafac=0;
    else
      gammafac=gamma_fourth(epstZ, tolerance);
    end
    tmp=0.4026*(tZ)^(  0.25);
    murow(i) = murow(i)  +  tmp*gammafac; 
    dmurow(i) = dmurow(i) + tmp*( ...
                    -epsf*exp(-epstZ)*(epstZ)^(-5/4) ...
                    +gammafac/(4*tZ))*Zi; 

    % 0.4026 = 0.306*3^(1/4) --> the factor 3 originates from 
    % 0.1342 = 0.306/3^(3/4) --> the factor 3 originates from 
    % the normalisation using the repation time td0 = 3*taue*Z^3

    % Scale the derivative
    dmurow(i)=-4*Z*(tZ)^(0.75)*dmurow(i)/3^0.25;
  end
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
