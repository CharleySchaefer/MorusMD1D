function calcChainExtension_B()
  close all;
  Z =10; % Number of entanglements
  Ns=11;
  We=0.7;
  Ntime=3000;
  

  R     =linspace(-Z*(Ns/(Ns-1))/2,Z*(Ns/(Ns-1))/2, Ns) % initial configuration
  dcorr=R(end)-R(1)-Z 
  ds0   = R(3)-R(2)
  dRds0 = R(2)-R(1)

  d2si=1.0/(pi^2*(ds0)^2)
  dt=0.4/d2si
  R(end)-R(1)-dcorr
  d2Rds2 = zeros(1,Ns);
  
  % Boundary condition
  time_current=0; lambda=1;
  tic
  for ti=1:Ntime
    if ti==Ntime/2
   %   We=0
    end
    %d2Rds2(1) =0; 
    %d2Rds2(Ns)=0;
  
    for i=2:Ns-1
      d2Rds2(i) = (d2si)*( R(i-1) - 2*R(i) + R(i+1) );
    end
    %            strain   convection
    R(2:Ns-1) = R(2:Ns-1) + dt*( d2Rds2(2:Ns-1) + (We/Z**2)*R(2:Ns-1) );
    R(1)  = R(2    ) - ds0;  % boundary condition
    R(end)= R(end-1) + ds0;  % boundary condition
    lambda=((R(end)   - R(1))-dcorr)/Z;
    R=R-mean(R)/2; % correct mean position
    time_current=time_current+dt;
    trow(ti) = time_current;
    lrow(ti) = lambda;
    sigma(ti)= get_stress(ds0, R);
  end
  toc
lambda
  lambda_local=zeros(1,Ns-2);
  for i=2:Ns-1
    lambda_local(i-1) = (R(i+1)-R(i-1))/(2*ds0);
  end


  figure
  subplot(2,2,1)
  loglog(trow/Z**2, lrow)
  xlabel('time/\tau_{R}')
  ylabel('\lambda')
  subplot(2,2,2)
  plot(linspace(-Z/2,Z/2, Ns), R, '.k', 'MarkerSize', 15); hold on
  plot(linspace(-Z/2,Z/2, Ns), sin( pi*sqrt((We/Z**2))*(  linspace(-Z/2,Z/2, Ns)  )  )...
                             /(     pi*sqrt((We/Z**2))*cos( pi*sqrt((We/Z**2))*(Z/2))  )  , 'r', 'LineWidth', 2)
  xlabel('s')
  ylabel('R')
  subplot(2,2,3)
  loglog((trow/Z**2), (sigma)); hold on
  xlabel('time/\tau_{R}')
  ylabel('\sigma')
  subplot(2,2,4)
  plot(linspace(-Z/2,Z/2,Ns-2)/Z, (lambda_local), '.k', 'MarkerSize', 15); hold on
  plot(linspace(-Z/2,Z/2,Ns)/Z, cos( pi*sqrt((We/Z**2))*(  linspace(-Z/2,Z/2, Ns)  )  )...
                             /(     cos( pi*sqrt((We/Z**2))*(Z/2))  )  , 'r', 'LineWidth', 2)
  xlabel('s')
  ylabel('\lambda')


  figure
  plot(linspace(-Z*(Ns/(Ns-1))/2,Z*(Ns/(Ns-1))/2, Ns), 2*ones(size(linspace(-Z*(Ns/(Ns-1))/2,Z*(Ns/(Ns-1))/2, Ns))), '.r', 'MarkerSize', 15); hold on;
  plot(R, ones(size(R)), '.k', 'MarkerSize', 15)
  plot(R(2:end-1), 1+1.0./lambda_local/Z, 'k', 'LineWidth', 3)
  plot(R(2:end-1), 1-1.0./lambda_local/Z, 'k', 'LineWidth', 3)
end

function sigma=get_stress(ds0, R)
  NR=length(R);
  dRds=zeros(1,NR);
  sigma=0;
  dRds(1)=0;dRds(end)=0;
  for i=2:NR-1
    dRds(i)=(R(i+1)-R(i-1))/(2*ds0);
  end
  for i=1:NR
    sigma = sigma + dRds(i)*dRds(i)*ds0;
  end

end

