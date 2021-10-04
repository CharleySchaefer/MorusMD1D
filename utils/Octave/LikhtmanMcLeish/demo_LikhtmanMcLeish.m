function demo_LikhtmanMcLeish()
  clc; close all;
  Z=300; Ntime=100;
  Ge=1.0;
  Nbeads=300;
  tolerance=1e-4;
  Cnu=1.0;
  trow=10.^linspace( 0, 10, Ntime  ); % units of taue
  taud0=3*Z^3;
 

  tic
  % function input in units of taud0
  mu=mu_LikhtmanMcLeish(    trow/taud0, Z, tolerance);
  R = R_LikhtmanMcLeish(Cnu*trow/taud0, Z, Nbeads);
  G = (4/5)*Ge*R.*mu;
  toc

  % Differentiate
  dRdt =zeros(1,Ntime-1);
  dmudt =zeros(1,Ntime-1);
  dtrow=zeros(1,Ntime-1); 
  trow_d=zeros(1,Ntime-1); 
  for i=1:Ntime-1
    trow_d(i) = (trow(i+1)+trow(i))/2;
    dtrow( i) = (trow(i+1)-trow(i));
    dRdt(  i) = (R(i+1)-R(i))/dtrow(i);
    dmudt(  i) = (mu(i+1)-mu(i))/dtrow(i);
  end
  dRdt_scaled =-4*Z*trow_d.^0.75.*dRdt;
  dmudt_scaled=-4*Z*trow_d.^0.75.*dmudt;

  ifp=fopen('demo_data.out', 'w');
  fprintf(ifp, '%12s %12s %12s %12s %12s\n', 'time/taue', 'mu(t)', 'dmu(t)', 'R(t)', 'dR(t)')
  for i=1:Ntime-1
    fprintf(ifp, '%12e %12e %12e %12e %12e\n', trow_d(i), mu(i), dmudt_scaled(i), R(i), dRdt_scaled(i));
  end
  fclose(ifp);

  figure
  subplot(1,2,1)
  loglog(trow(dmudt_scaled>1e-8), dmudt_scaled(dmudt_scaled>1e-8)); hold on;
  loglog(trow(dRdt_scaled>1e-8), dRdt_scaled(dRdt_scaled>1e-8) );
  xlim( [1, 1e10]); ylim([0.1, 4]); %([0,0.1, 1e10,  4])
  xlabel('t/\tau_e'); ylabel('dR, d\mu')
  legend('\mu(t)', 'R(t)', 'Location', 'southwest')
  subplot(1,2,2)
  loglog(trow(G>1e-8), G(G>1e-8));
  xlabel('t/\tau_e'); ylabel('Relaxation Modulus')
end
