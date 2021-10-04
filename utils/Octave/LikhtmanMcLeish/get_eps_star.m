% Solve Eq.(14) in Likhtman-McLeish
%  given in units of taud0*Z, with taud0=3*taue*Z^3;
%  the factor 0.403 = 3^(1/4)*0.306, with 0.306 the value in the paper
function eps_star=get_eps_star(Z,Gf)
  p_star=sqrt(Z*0.1);
  dsum=0.0; p=1;
  while p<=p_star
    dsum=dsum+1.0/(p*p);
    p=p+2;
  end
  eps_star=(4*0.403/(1-dsum*8*Gf/pi^2))^4;
end
