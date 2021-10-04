% Verification: checkFig3.m
function Gf=get_renormalisation_elastic_modulus(Z)
  sqZi=1.0/sqrt(Z);
  Gf=1.0+(-1.69 + (2.0 - 1.24*sqZi)*sqZi)*sqZi;
end
