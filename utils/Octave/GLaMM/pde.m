% GLaMM
function pde()
  clc; close all;
  % p=1,2,...,N+1
  % q=1,2,...,N+1

  %=====================================
  % USER INPUT
  flow_mode='extensional'
  strain_rate=0.01;
  N_num=20;
  a=1.0; % tube diameter
  ai=1.0/a;
  Rs=2.0;
  tolY=4e-5;
  NtimeMax=400;

  tau_e=1.0;
  Z=3; % number of entanglements

  % Adaptive step size / time stepper
  delta_t=0.8*tau_e*(Z*1.0/N_num)^2; % first time step
  delta_t0=delta_t;
  tol=0.2; % tolerance


  % END USER INPUT
  %=====================================


  
  %=====================================
  % INITIALISE
  Kappa=zeros(3); % flow tensor
  if(flow_mode=='extensional')
    Kappa(1,1) = -0.5*strain_rate;
    Kappa(2,2) = -0.5*strain_rate;
    Kappa(3,3) =  1.0*strain_rate;
  end

  dsi=N_num/(2*Z);
  d2si=(N_num/Z)^2;
  
  feq=zeros(N_num+1, N_num+1, 3, 3);
  for p=1:N_num+1
    for q=1:N_num+1
      if abs(p-q)<dsi
        feq(p,q,1,1)=1.0/3;
        feq(p,q,2,2)=1.0/3;
        feq(p,q,3,3)=1.0/3;
       end
    end
  end

  f=feq;
  t=0;
  dfdt=zeros(N_num+1, N_num+1, 3, 3);

%  d2mat=zeros(size(dfdt));
  prefac_reptation  = d2si*(1.0/(3*pi^2*Z*tau_e));
  prefac_retraction = Rs/(2*pi^2*tau_e);
 %dfdt=zeros(N_num+1, N_num+1, dsi); 

trow=[];
errYrow=[];
sigma_arr=[];
Zstarrow=[];
measure1=[];
measure2=[];
measure3=[];
sigma_arr=zeros(NtimeMax, 3,3);
dtrow=[];
diagonal=zeros(N_num+1,1);

  % INITIALISATION DONE
  %=====================================
  Zstar = fnc_Zstar(f, ai, Z,N_num); 
  sigma=zeros(3,3);
  for i=1:N_num+1
  for m=1:3
  for n=1:3
    sigma(m,n) = sigma(m,n)+f(i,i,m,n);
  end
  end
  end
  sigma=sigma*12/(5*(N_num+1));

  fprintf("%12s %12s %12s %12s %12s %12s\n", "#time", "dt", "Zstar", "sigma_xx", "sigma_yy", "sigma_zz");
  fprintf("%12e %12e %12e %12e %12e %12e\n", 0.0, 0.0, Zstar, sigma(1,1), sigma(2,2), sigma(3,3));

Ntime=NtimeMax;
for ti=1:NtimeMax

  % Convection term
  if(flow_mode=='extensional')
    for p=1:3
      for q=1:3
        dfdt(:,:, p,q) = ( Kappa(p,p) + Kappa(q,q) )*f(:,:,p,q);
      end
    end
  end
  %fprintf("testA: %e\n", sum(dfdt(:)))
  %measure2=[measure2,sum(dfdt(:))];

  % Reptation term
  Zstar = fnc_Zstar(f, ai, Z,N_num); 
  prefact_rep_tmp=prefac_reptation*(Z/Zstar)^2;
  for p=2:N_num
    for q=2:N_num
      dfdt(p,q,:,:) = dfdt(p,q,:,:) + ...
       prefact_rep_tmp*( f(p+1,q+1, :,:)+f(p-1,q-1, :,:)-2*f(p,q, :,:) );
    end
  end
  %fprintf("testB: %e\n", sum(dfdt(:)))

  % Retraction term
  %Tr_f=trace4(f);
  for q=1:N_num+1
    diagonal(q)=log( f(q,q,1,1)+f(q,q,2,2)+f(q,q,3,3) );
  end 

  ddiag=zeros(N_num+1,1);
  for q=2:N_num
    ddiag(q)=dsi*(diagonal(q+1)-diagonal(q-1));
  end 

  retraction=prefac_retraction* ...
    (
      dmat_ds(  prod41B(f, ddiag), N_num, dsi ) + ...
      dmat_dsB( prod41A(f, ddiag), N_num, dsi )...
    );


  dfdt = dfdt + retraction;
  %fprintf("testC: %e\n", sum(dfdt(:)))
  %measure3=[measure3,sum(retraction(:))];

  %measure1=[measure1,sum(dfdt(:))];
  % Update 


  t=t + delta_t; % update time


  for m=1:3 % boundary condition
    for n=1:3
      dfdt(1,1,m,n)=0;
      dfdt(N_num+1,N_num+1,m,n)=0;
    end
  end
  deltaf=delta_t*dfdt;
  f=f + deltaf;

  % update time step
  err=max(abs(deltaf(:)./f(:)));
  delta_t=delta_t*(tol/(0.02*tol+err))^0.25;
  delta_t=max(delta_t, delta_t0);
  % factor 0.02*tol prevents too large jumps in time step
  

  sigma=zeros(3,3);
  for i=1:N_num+1
  for m=1:3
  for n=1:3
    sigma(m,n) = sigma(m,n)+f(i,i,m,n);
  end
  end
  end
  sigma=sigma*12/(5*(N_num+1));




  trow=[trow, t];
  dtrow=[dtrow, delta_t];
  sigma_arr(ti,:,:)= sigma;
  Zstarrow=[Zstarrow, Zstar];

  if ti>10
    tmp2 = sigma_arr(ti,:,:);
    tmp  = sigma_arr(ti-10,:,:);
    errY = max(abs((tmp(:)-tmp2(:))./(1e-12+tmp(:)) )) *10*delta_t0/(trow(ti)-trow(ti-10)) ;
    if(errY<tolY)
      Ntime=ti-1;
      break;
    end
  else
    errY=tolY;
  end
  errYrow=[errYrow, errY];
  fprintf("%12e %12e %12e %12e %12e %12e\n", trow(ti), dtrow(ti), Zstarrow(ti), sigma_arr(ti, 1,1), sigma_arr(ti, 2,2), sigma_arr(ti, 3,3));

end

delta_t

sigma_arr(Ntime,:,:)
  figure
  subplot(2,2,1);
 plot(trow(1:Ntime), log10(sigma_arr((1:Ntime),1,1))); hold on
 plot(trow(1:Ntime), log10(sigma_arr((1:Ntime),2,2))); hold on
 plot(trow(1:Ntime), log10(sigma_arr((1:Ntime),3,3))); hold on
 plot(trow(1:Ntime), log10(sigma_arr((1:Ntime),1,2))); hold on
 plot(trow(1:Ntime), log10(sigma_arr((1:Ntime),1,3))); hold on
 plot(trow(1:Ntime), log10(sigma_arr((1:Ntime),2,3)));
  ylabel('sigma')
  xlabel('time')


  subplot(2,2,2);
 plot(trow(1:Ntime), Zstarrow(1:Ntime));
  ylabel('Z')
  xlabel('time')
%  subplot(2,2,3);
% plot(trow, measure1); hold on
% plot(trow, measure2); hold on
% plot(trow, measure3); hold on
% legend('all', 'flow', 'retraction')
  subplot(2,2,3);
  plot(trow(1:Ntime), log10(dtrow(1:Ntime))); hold on
  plot(trow(1:Ntime), log10(delta_t0*ones(size(trow(1:Ntime)))) );
  ylabel('time stepsize')
  xlabel('time')

  subplot(2,2,4);
  plot(trow(1:Ntime), 1:Ntime); hold on
  ylabel('number of time steps')
  xlabel('time')

  figure
  plot(trow(1:Ntime), log10( errYrow(1:Ntime))); hold on;
  plot(trow(1:Ntime), log10( tolY*ones(size(trow(1:Ntime))))); hold on;



end

function out=prod41A(mat4, mat1)
  N=length(mat1);
  out=zeros(N,N,3,3);
for p=1:N
for q=1:N
  for m=1:3
    for n=1:3
      out(p,q,m,n) = mat4(p,q,m,n).*mat1(p);
    end
  end
end
end
end

function out=prod41B(mat4, mat1)
  N=length(mat1);
  out=zeros(N,N,3,3);
for p=1:N
for q=1:N
  for m=1:3
    for n=1:3
      out(p,q,m,n) = mat4(p,q,m,n).*mat1(q);
    end
  end
end
end
end

function out=prod42(N_num, mat4, mat2)
  out=zeros(N_num+1,N_num+1,3,3);
  for m=1:3
    for n=1:3
      out(:,:,m,n) = mat4(:,:,m,n).*mat2;
    end
  end
end

function mat2=trace4(mat4)
  mat2=mat4(:,:,1,1)+mat4(:,:,2,2)+mat4(:,:,3,3);
end

function Zstar = fnc_Zstar(f, ai, Z,N_num)
  Zstar = 0;
  for p=1:(N_num+1)
    Zstar = Zstar + sqrt( f(p,p,1,1)+f(p,p,2,2)+f(p,p,3,3) );
  end
  Zstar=Zstar*ai*Z/(N_num+1);
end

function dmat = dmat_ds(  mat, N_num, dsi)
  dmat=zeros(size(mat));
  for p=2:N_num
    dmat(p,:,:,:)=dsi*(   mat(p+1,:,:,:)-mat(p-1,:,:,:) );
  end
end
function dmat = dmat_dsB( mat, N_num, dsi)
  dmat=zeros(size(mat));
  for q=2:N_num
    dmat(:,q,:,:,:)=dsi*(   mat(:,q+1,:,:)-mat(:,q-1,:,:) );
  end
end

% Eq. (42b) in Graham et al. (2003)
function d2mat = d2mat_ds2(mat, d2si)
  d2mat=zeros(size(mat));
  for p=2:N
    d2mat(p,:,:,:)=d2si*( mat(p+1,:,:,:)+mat(p-1,:,:,:) -2*mat(p,:,:,:) );
  end
end

% Eq. (43a) in Graham et al. (2003)
function dmat = dmat_clf(mat, dsi)
  dmat=zeros(size(mat));
  for p=2:N
    for q=2:N
      dmat(p,q,:,:)=dsi*(mat(p+1,q+1,:,:)-mat(p-1,q-1,:,:));
    end
  end
end

% Eq. (43b) in Graham et al. (2003)
function d2mat = d2mat_clf(mat, N_num, d2si)
  d2mat=zeros(size(mat));

  for p=2:N_num
    for q=2:N_num
      d2mat(p,q,:,:)=d2si*( mat(p+1,q+1, :,:)+mat(p-1,q-1, :,:)-2*mat(p,q, :,:) );
    end
  end
end
