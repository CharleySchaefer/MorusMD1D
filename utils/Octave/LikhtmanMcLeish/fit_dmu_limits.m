% Extract parameters from simulation
% Split data in early stage [taue - tauR]     --> Cmu
% and late stage            [tauL - tau_end ] --> Gf and Tf
function [Cmu, Gf, Tf]=fit_dmu_limits(fname, ZE, run_in_matlab)
  clc; close all;
  %=========================================
  % INPUT
  % columns 'time' and 'dmuII' will be used for curve fitting
  %fname='mu.out'; % columns: time, muI, dmuI, muII, dmuII 
  %run_in_matlab=0; 
  %ZE=10;
  %=========================================


  %=========================================
  % IMPORT DATA 
  delim=' ';

  data=importdata(fname);
  Nheader=1;
  c=data(Nheader,1){1}(1);
  while c<'0' || c>'9'
    Nheader=Nheader+1;
    c=data(Nheader,1){1}(1);
  end
  Nheader=Nheader-1;
  data=importdata(fname, delim, Nheader);

  data=data.data;
  tdata=data(:, 1);
  mudata=data(:,5);
  Ndata=length(tdata);
  tauR=ZE*ZE; % Separation between early and late times
  %=========================================


  %=========================================
  % FIT EARLY-TIME PLATEAU
  indR=0;
  ind0=1; tt=tdata(ind0); countI=0; countII=0;
  dmuI=0; Cmu=0;
  while tt<tauR
    if indR==0 && tt>=1
      indR=ind0;
    end
    tmp=mudata(ind0);
    if tmp>0 && tt>=1
      Cmu=Cmu+tmp; countII=countII+1;
    end
    ind0=ind0+1; tt=tdata(ind0);
  end
  Cmu=Cmu/countII;
  %=========================================


  %=========================================
  % FIT LATE-TIME PLATEAU
  IncludeOptimisationPkg(run_in_matlab); % include optimisation package

  % Get data range tauL-tau_end
  indLateU=Ndata;
  while isnan( tdata(indLateU) ) || isinf( tdata(indLateU) )
    indLateU=indLateU-1;
  end
  indLateL=indLateU-1;

  while mudata(indLateL)>mudata(indLateL+1) || mudata(indLateL)<Cmu
    indLateL=indLateL-1;
  end
  while mudata(indLateL)<mudata(indLateL+1) 
    indLateL=indLateL-1;
  end
  t_late =tdata(indLateL:indLateU);
  dmu_late=data(indLateL:indLateU, 5);

  [param, param_std, Rsquared] = FitMuLate(t_late, dmu_late, ZE);
  Gf=param(1);
  Tf=param(2);
  %=========================================

  %=========================================


  %=========================================
  % REPORT  
  if 0
  fprintf('Early-stage fit:\n');
  fprintf(' Cmu=%f\n', Cmu);
  fprintf('Late-stage fit (R^2 = %f).\n', Rsquared);
  fprintf('  Gf: %e +/- %e\n', param(1), param_std(1));
  fprintf('  Tf: %e +/- %e\n', param(2), param_std(2));
  
  mu_late=cost_mu_late( param, ZE, t_late, 0, 1e-6) ; %mu_late.*ss;
  yfit=cost_dmu(paramF, ZE, tfull/(3*ZE*ZE*ZE), 0,  1e-5);

  figure
  plot(log10(tdata), data(:, 3), '.r', 'MarkerSize', 15) ; hold on
  plot(log10(tdata), data(:, 5), '.b', 'MarkerSize', 15) ; hold on
  %plot( log10([tdata(1), tauR]), dmuI*[1 1], 'r', 'LineWidth', 2)
  plot( log10([tdata(1), tauR]), Cmu*[1 1], 'k', 'LineWidth', 2)
  plot( log10([tdata(1), tauR]), Cmu*[1 1], 'k', 'LineWidth', 2)
  plot( log10([t_late]), mu_late, 'k', 'LineWidth', 2)
  plot( log10(tfull), yfit, 'g', 'LineWidth', 2)
  axis( [ log10(tdata(1)), log10(6*3*ZE*tauR), -2, 5] )
  xlabel('log_{10}(t/\tau_e)')
  ylabel('d\mu')
  end
end

function cost=cost_mu_late( param, ZE, tdata, dmudata, tol)
  Gf=param(1);
  Tf=param(2);
  td0=3*ZE^3; % reptation time
  dmu_fnc=32*Gf*tdata.^(0.75)/(3*Tf*(pi*ZE)^2);
  ss=zeros(size(tdata)); err=2*tol;
  p=1;
  
  while err>tol
    dss=exp(-p*p*tdata/(Tf*td0) );
    ss=ss+dss;
    err=dss/ss;
    p=p+2;
  end
  dmu_fnc=dmu_fnc.*ss;
  cost=dmu_fnc-dmudata;
end

function IncludeOptimisationPkg(run_in_matlab)
  %====================================================
  % INCLUDE OPTIMISATION PACKAGES
  %----------------------------------------------------
  if run_in_matlab % RUN IN MATLAB
    % Check licences
    if ~license('test','optimization_toolbox')
      error('ERROR: this code uses MATLAB''s optimization toolbox!')
      return
    end
    if ~license('test','statistics_toolbox')
      error('ERROR: this code uses MATLAB''s statistics toolbox!')
      return
    end
  else % RUN IN OCTAVE
    % install octave packages: run "pkg install -forge struct optim stk statistics"
    %                          in command window
	try   pkg load struct
	catch
	  printf("To install struct packakge, run pkg install -forge struct."); exit
	end
	try   pkg load optim
	catch
	  printf("To install struct packakge, run pkg install -forge optim."); exit
	end
%    pkg load statistics  % lhsdesign (not yet implemented)
  end
  % OPTIMISATION PACKAGES INCLUDED
  %====================================================
end

%====================================================
% FITTING ALGORITHM
%   Input
%     tdat: array with times 
%     pdat: array with fractions of cells with an aggresome at time tdat 
%   Output
%     param:     [n,J] 
%     param_std: uncertainty in [n, J]
%     Rsquared:  Fit quality
function [param, param_std, Rsquared] = FitMuLate(tdata, dmudata, ZE)
  tol=1e-6;
  %----------------------------------------------------
  % FIT
  options = optimset(...   
            'MaxIter',1000,...
            'Display','off',...
            'MaxFunEvals',100000,...
            'TolX',1e-10,...
            'TolFun',1e-10);%,...
%             'Algorithm','levenberg-marquardt'); % 'trust-region-reflective'

%	param0=[ones(size(xdat)), xdat]\ydat ;% linear regression
  param0=[1,1];

[param, resnorm, residual, qqq1, qqq2, qqq3, jacobian]=...
       lsqnonlin(@(param)cost_mu_late( param, ZE, tdata, dmudata, tol), param0, [0 0], [1 1], options);

  %----------------------------------------------------
  % ERROR ANALYSIS

  chisquare= resnorm; 
  Npar = 2;                              % Number of fit parameters
  dof  = length(tdata) - Npar;            % Degrees of freedom
p=chi2cdf(chisquare, dof);
  C          = inv(full(jacobian)'*full(jacobian));
  covarpar   = chisquare/dof*C;            % Variance-covariance matrix
  param_std  = sqrt(diag(covarpar))';    % Standard deviation (STD)

  corr_mat = C./sqrt(diag(C)*diag(C)');  % Correlation matrix 

  Rsquared=1-resnorm/((dmudata-mean(dmudata))'*(dmudata-mean(dmudata)));
end
