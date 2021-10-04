% Verification: checkFig3.m
function Tf=get_renormalisation_reptation_time(Z)
  sqZi=1.0/sqrt(Z);
  Tf=(1+sqZi*(-2*1.69+sqZi*(4.17-1.55*sqZi)));
end
