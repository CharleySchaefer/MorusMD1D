% MODEL REF: Leibler,Rubinstein & Colby, Macromolecules 24, 4701--4707 (1991)
% Implementation: C. Schaefer 2020
% Output: diffusivity in units of (non-sticky) Rouse diffusivity (Drouse = kT/(zeta*N))
function D=StickyDiffusionCoefficient(p, ZE, ZS, tauS)
  tauR=ZE*ZE; % Rouse time of the chain, in units of the Rouse time of an entanglement strand

  % Get kmax
  if p>0
    MAXFAC=( tauS*(ZS*1.0/ZE)^2 )*(1.0-p)/p;
    arr=roots([1, 2, 1, -MAXFAC]);
    kmax=arr(3);
  elseif p==0
    kmax=inf;
  end
  
  % Sum over k's
  Dksum=0; Dkend=0;
  for k=1:ZS-1
    if k<=kmax
      Dk_fac=(k+1)^3;
    else
      Dk_fac=(kmax+1)^3*sqrt(kmax*1.0/k);
    end
	
    fac=k*(1-p)^(k-1)*Dk_fac;
    if k<= ZS-2
      Dksum=Dksum+(ZS-k-1)*fac;
    end
    Dkend=Dkend+fac;
  end
  fac=p^2/(tauS*(ZS+1)^3); % NOTE: LRC paper includes an erroneous factor 1/4 
  Dksum=  p*Dksum*fac;
  Dkend=4.5*Dkend*fac;
	
  if tauS<=tauR
    DS_fac=(kmax+1)*1.0/(ZS+1)*sqrt(ZS*1.0/kmax);
  else
    DS_fac=1.0;
  end
  DS=(1.0/tauR)*(1-p)^ZS*DS_fac;

  D=( Dksum+Dkend+DS)*tauR;
end
