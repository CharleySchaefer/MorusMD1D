%--------------------------------------------------------
% INPUT: 
%   trow:      array  - Time in units of tau_e 
%   Z:         double - Number of entanglements per chain
%   Cnu:       double - strength of countour-length fluctuations 
%   tolerance: double - numerical setting - tested for tolerance=1e-3
% OUPUT:
%   Gt:        array - Relaxation modulus (in units of Ge)
function Gt=Gt_LikhtmanMcLeish(trow, Z, Cnu, tolerance)
  tau_d0 = 3*Z^3; % Reptation time
  murow  = mu_LikhtmanMcLeish(trow/(3*Z^3), Z, tolerance);
  Rrow   = R_LikhtmanMcLeish( trow, Z, Cnu);
  Gt     = 0.8*murow.*Rrow;
end
