%-------------------------------------------------------
% PURPOSE
% Calculate Likhtman-McLeish's mu(t).
%   Coded by C. Schaefer (2020)
%   Verified by comparison of Ze=100 output to the curve 
%   in Fig 4 in Likhtman-McLeish paper
%-------------------------------------------------------
function get_mu_LMtheory(varargin)
%-------------------------------------------------------
% USER INPUT
  Narg=length(varargin);

  if Narg<4
    printf('at least four arguments expected.\n');
    printf('arguments: Ze, tL, tU, tbins, (optional)tol\n')

    break;
  end

                                     % settings 
  Ze=(varargin{1});           % 1. Number of entanglements
  tL=(varargin{2})/(3*Ze^3);  % 2. Lower bound time interval 
  tU=(varargin{3})/(3*Ze^3);  % 3. Upper bound time interval
  nbins=(varargin{4});        % 4. Number of bins in time interval
  if Narg>4
    tolerance=(varargin{5});  % 5. Numberical tolerance
  else
    tolerance=1e-4;
  end
%-------------------------------------------------------

%-------------------------------------------------------
% INCLUDE LIKHTMAN-MCLEISH MODULE
mfile_callers=dbstack;
fname=mfile_callers.file;
while( fname(end) ~= '/' )
  fname(end)=[];
end
dname=strcat(fname, 'LikhtmanMcLeish');
addpath(dname);
%-------------------------------------------------------

%-------------------------------------------------------
% CALCULATE
  % Time range (units of reptation time)
trow=10.^linspace(log10(tL), log10(tU), nbins  ) ;

% Calculate mu(t) on time range
[mu dmu]=mu_LikhtmanMcLeish(trow, Ze, tolerance);

%-------------------------------------------------------

%-------------------------------------------------------
% EXPORT
% Export t and mu(t) columns
for i=1:length(mu)
%  if i==1
    mu0=mu(i); ttime0=trow(i);
    fprintf('%12.3e %12.3e %12.3e\n', (3*Ze^3)*ttime0, mu0, dmu(i));
%  else % TO COMPARE dmu(i) TO NUMERICAL VALUE
%    mu1=mu(i);
%    ttime=trow(i); 
%    dmu_num=-4*Ze*((3*Ze^3)*ttime)^0.75*(mu1-mu0)/((3*Ze^3)*(ttime-ttime0));
%    fprintf('%12.3e %12.3e %12.3e %12.3e\n', (3*Ze^3)*trow(i), mu(i), dmu_num, dmu(i));
%    mu0=mu1; ttime0=ttime;
%  end
end
  
end
