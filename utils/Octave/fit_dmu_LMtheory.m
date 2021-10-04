%-------------------------------------------------------
% PURPOSE
% Calculate Likhtman-McLeish's mu(t).
%   Coded by C. Schaefer (2020)
%   Verified by comparison of Ze=100 output to the curve 
%   in Fig 4 in Likhtman-McLeish paper
%-------------------------------------------------------
function fit_dmu_LMtheory(varargin)

%-------------------------------------------------------
% USER INPUT
  nargin;
  Narg=length(varargin);

  if Narg<3
    printf('Three arguments expected.\n');
    printf('arguments: fmuname, ZE, run_in_matlab\n')

    break;
  end
                                     % settings 
  fmuname=varargin{1};
  ZE=(varargin{2});           % 1. Number of entanglements
  run_in_matlab=(varargin{3}); % 2. Lower bound time interval 
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
% FIT
[Cmu, Gf, Tf]=fit_dmu_limits(fmuname, ZE, run_in_matlab);
fprintf('#Cmu Gf Tf\n')
fprintf('%f %f %f\n', Cmu, Gf, Tf);
%-------------------------------------------------------
end

  
