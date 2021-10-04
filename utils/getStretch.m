 function getStretch(fname, ZE, Nbins)
  clc; close all;

  [lambda, P]=getStretchDistribution(fname, ZE, Nbins);

  for i=1:Nbins
    fprintf('%12e %12e\n', lambda(i), P(i));
  end
end

function [lambda, P]=getStretchDistribution(fname, ZE, Nbins)

   data=importdata(fname, ' ', 1);

  data=data.data;
  L=abs(data(:,3)-data(:,2))/ZE;


  lL=log10(L);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % BIN ON LOG SCALE - HISTCOUNTS
  % [N,edges] = histcounts(lL, Nbins); % Not yet implemented in octave
  Lmin=min(lL);
  Lmax=max(lL);
  edges=linspace(Lmin-1, Lmax+1, Nbins+1); dL=edges(2)-edges(1);
  N=zeros(1,Nbins);
  for i=1:length(L)
    ind=round( (lL(i)-(Lmin-1) )/dL ) + 1;
    N(ind)=N(ind)+1;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % NORMALISE ON LINSCALE
  lambda=zeros(1,length(N));
  norm_fac=sum(N);
  for i=1:length(N)
    lambda(i)=10^(edges(i)+0.5*dL);
    P(i)=N(i)*1.0/( 10^(edges(i+1)) - 10^(edges(i)) )/norm_fac;
  end



end
