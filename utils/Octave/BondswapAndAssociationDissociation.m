% Model influence of bondswapping on sticker opening rate
function BondswapAndAssociationDissociation()
close all; clc;
  X=10.^(linspace(-3,3,20));
  p=pfnc(X);

  Yrow=[0.1 0.2 0.5 1 2 5 10 20 50 100];

  figure 
  for i=1:length(Yrow)
	  Y=Yrow(i);
  Ropen  =1+Y.*(1-p);
  Rclose =X.*(1-p) + 0.5*Y*p;

  if Y==1
  plot(log10(X), log10(Ropen), 'k'); hold on;
  plot(log10(X), log10(Rclose), 'k');
  else
  plot(log10(X), log10(Ropen), 'r'); hold on;
  plot(log10(X), log10(Rclose), 'b');
end
  end
  legend('Rate of sticker opening', 'Rate of sticker closing')
  xlabel('fraction of closed stickers');
  ylabel('log Rate (in units of bare dissociation rate)')'
 % plot(log10(X), pfnc(X))
  

end

function p=pfnc(X)
  p=1-(-1+sqrt(1+8*X))./(4*X);
end
