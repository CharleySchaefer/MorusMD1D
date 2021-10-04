function AnalyseRe(pth,Nbins)
  close all; clc; fclose all;
  %=======================
  % USER INPUT
%  Nbins=20;
%  pth='.';
  %=======================

  f_in = sprintf('%s/timeprogress.out',pth);
  delim=' ';
  Nheader=1;
  col1=10;
  binfac=2;
  data=importdata(f_in, delim, Nheader);
  
  data=data.data;
  [Ntime, col2]= size(data);
  
  Ncol=col2-col1+1;

  t_dat =data(:,1);
  Re_dat=data(:,col1:col2);


  Re_min= min(Re_dat(:));
  Re_max= max(Re_dat(:));

  % Binning on log scale
  lRe_min=log10(Re_min/binfac);
  lRe_max=log10(binfac*Re_max);
  edges=linspace(lRe_min, lRe_max, Nbins+1);
  dlRe=edges(2)-edges(1);

  Re=zeros(1,Nbins);
  P =zeros(Ntime,Nbins);
  for i=1:Ntime
    bins=zeros(1,Nbins);
    lRe_row=log10(Re_dat(i,:));
    for j=1:Ncol
      ind= 1+round( (lRe_row(j)-lRe_min)/dlRe );
      bins(ind)=bins(ind)+1;
    end

    % Normalise 
    for j=1:Nbins
      Re(j)=10^(edges(j)+0.5*dlRe);
      P(i,j)=bins(j)*1.0/( 10^(edges(j+1)) - 10^(edges(j)) )/Ncol;
    end
  %  if i==1
  %    plot(Re, P, 'k'); hold on
  %  else
  %    plot(Re, P); hold on
  %  end
  end
  %=======================
  % Export
  f_out= sprintf('%s/stretch_distribution_transient.out',pth);
  ifp=fopen(f_out, 'w');
  for j=1:Nbins
    fprintf(ifp, '%12e ', Re(j));
    for i=1:Ntime
      fprintf(ifp, '%12e ', P(i,j));
    end
    fprintf(ifp, '\n');
  end
  fclose(ifp);
  Pmean=zeros(1,Nbins);
  Pvar =zeros(1,Nbins);
  for j=1:Nbins
    Pmean(j)=mean(P(:,j)); % Time average
  end
  for j=1:Nbins
    Pvar(j)=mean( (P(:,j)- Pmean(j)).^2 );
  end

  f_out= sprintf('%s/stretch_distribution_steady.out',pth);
  ifp=fopen(f_out, 'w');
  
  fprintf(ifp, "%12s %12s %12s\n", 'Re', 'Pmean', 'Psem' );
  for j=1:Nbins
    fprintf(ifp, "%12e %12e %12e\n", Re(j), Pmean(j), sqrt(Pvar(j)/(Ncol-1) ) );
  end

  fclose(ifp);
  
end
