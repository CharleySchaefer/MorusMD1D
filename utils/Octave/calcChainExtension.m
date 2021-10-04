function calcChainExtension()
  clc; close all;
  Ze=10;    % Number of entanglements per chain

  Ns=100; % Discretisation of chain coordinate
  
  figure

  We=0.20;  % Weissenberg number
  for We=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]
  % Sticker distribution
  stickers=zeros(1,Ns);
  ns=5;
  while ns<=Ns
    stickers(ns)=1;
    ns=ns+10;
  end
  Zs=sum(stickers);

  % Chain coordinates on chain
  s_arr  = linspace(0, 1, Ns);

  % Chain coordinates in real space
  R_arr0=linspace(-Ze/2, Ze/2, Ns); % Quiescent 
  R_arr =             Ze*sin( [2*s_arr-1]*pi*sqrt(We)/2 )...
           /(pi*sqrt(We)*cos(             pi*sqrt(We)/2 ) );

  % Local chain stretch (dR/ds)
  lambda_arr = cos( [2*s_arr-1]*pi*sqrt(We)/2)/...
               cos(             pi*sqrt(We)/2);

  % Sticker density
  rho_s=zeros(1,Ns);
  S_arr=R_arr(find(stickers==1)); % Sticker coordinates
  L_arr=lambda_arr(find(stickers==1)); % Local stretch at sticker
  for i=1:Ns
    for ns=1:Zs
      varL=1.0/L_arr(ns)/3;
      rho_s(i)=rho_s(i)+exp(-(R_arr(i)-S_arr(ns))^2/(2*varL) )/sqrt(2*pi*varL);
    end
  end

  x=linspace(R_arr(1), R_arr(end), Ns);
  y0=interp1(R_arr, rho_s, x); 
  S0=round(interp1(R_arr, (stickers), x)); 
  yFT=(fft(y0));
  yFT=(fft(S0));
  y  = fftshift(abs( ifft( yFT.*S0 ) )); %/(2*Zs);

    subplot(2,2,1)
  plot(s_arr, R_arr  , 'r', 'LineWidth', 2); hold on;
  xlabel('chain coordinate, s')
  ylabel('spatial coordinate, R')

    subplot(2,2,2)
  plot(R_arr/max(R_arr), We*ones(length(R_arr)), '.k', 'MarkerSize', 5); hold on
  plot(R_arr(find(stickers==1))/max(R_arr), We*find(stickers==1)./find(stickers==1), '.r', 'MarkerSize', 10);; hold on
  plot(R_arr/max(R_arr), We+rho_s*0.05, 'g', 'LineWidth', 3) ; hold on
  %[AX,H1,H2] = plotyy(R_arr/max(R_arr), We+rho_s*0.05, R_arr/max(R_arr), We+rho_s*0.05); hold on;
 % plotyy(R_arr/max(R_arr), We+rho_s*0.05, R_arr/max(R_arr), We+rho_s*0.05); hold on;
  %[AX,H1,H2] = plotyy(x,y1,x,y2,'plot');
  xlabel('chain coordinate, s')%
%ylabel(AX(1),'Wi')
%ylabel(AX(2),'Sticker density')

    subplot(2,2,3)
  plot(s_arr, lambda_arr  , 'r', 'LineWidth', 2); hold on;
  xlabel('chain coordinate, s')
  ylabel('local extension, \lambda')

    subplot(2,2,4)
    plot(x,S0); hold on;
    plot(x,y)

%  ylabel('Wi')
  %plot(R_arr0, 0.1+ones(length(R_arr0)), '.k', 'MarkerSize', 10); hold on
  %plot(R_arr0(find(stickers==1)), 0.1+find(stickers==1)./find(stickers==1), '.r', 'MarkerSize', 15);
  end

end
