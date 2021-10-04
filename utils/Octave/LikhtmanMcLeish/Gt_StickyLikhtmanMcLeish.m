%--------------------------------------------------------
% INPUT: 
%   trow:      array  - Time in units of tau_s 
%   Z:         double - Number of entanglements per chain
%   Cnu:       double - strength of countour-length fluctuations 
%   tolerance: double - numerical setting - tested for tolerance=1e-3
% OUPUT:
%   Gt:        array - Relaxation modulus (in units of Ge)
function Gt=Gt_StickyLikhtmanMcLeish(trow, Ze, Zs, Cnu, tolerance)
  tau_d0=3*Zs^2*Ze; % Reptation time
  murow=mu_LikhtmanMcLeish(trow/tau_d0, Zs, tolerance);
  Rrow =R_LikhtmanMcLeish( trow, Zs, Cnu);
  Gt =0.8*murow.*Rrow;
end
