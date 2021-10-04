function average_timeprogress(pth)
%  close all; clc;

  Nheader=1; % Number of header lines in timeprogress files
  ycol=4;    % RCM - center of mass
%  if  nargin==0
%    printf("Usage: ./average_timeprogress <dir> <ycol to average>\n")
%    exit 
%  end

  % Get sudirectories with time progress files
%  arg_list = argv();
%  pth=arg_list{1};
  if pth(end)=='/'
    pth(end)=[];
  end

  subdirs=dir(sprintf('%s', pth));
  Ncontents=length(subdirs);
  for i=1:Ncontents
    if ( subdirs(Ncontents+1-i).isdir==1) && strncmp(subdirs(Ncontents+1-i).name, "seed", 4)
; % do nothing
    else % Remove entry
      subdirs(Ncontents+1-i)=[];
    end
  end
  Nseeds=length(subdirs);


  % AVERAGE
  for i=1:Nseeds
    rseed=i; %rseed_row(i);
    fname=sprintf('%s/%s/timeprogress.out', pth, subdirs(i).name);

    data=importdata(fname, ' ', Nheader); data=data.data;
    if i==1
      xdata=data(:,1);
      Ldata=(data(:,3)-data(:,2)); % Chain length
      ydata=data(:,ycol).^2;       % Rcm
      stress=(data(:,5));          % Stress

    else
    xdata=xdata+data(:,1);
    Ldata=Ldata+(data(:,3)-data(:,2)); % Chain length
    ydata=ydata+data(:,ycol).^2;       % Rcm
    stress=stress+(data(:,5));         % Stress
    end
  end
  xdata=xdata/Nseeds;
  ydata=ydata/Nseeds;
  ymean=sqrt(ydata);
  Ldata=Ldata/Nseeds;
  stress=stress/Nseeds;

  % STD
  for i=1:Nseeds
    rseed=i; %rseed_row(i);
    fname=sprintf('%s/%s/timeprogress.out', pth, subdirs(i).name);

    data=importdata(fname, ' ', Nheader); data=data.data;
    if i==1
      ystd=     (abs( data(:,ycol) ) - ymean).^2; % Rcm
      Lstd=     (  Ldata - (data(:,3)-data(:,2))    ).^2; % chain length
      stressstd=     (  stress - (data(:,5))    ).^2; % chain length
    else
      ystd=ystd+(abs( data(:,ycol) ) - ymean).^2; ; % Rcm
      Lstd=Lstd+(  Ldata - (data(:,3)-data(:,2))    ).^2; % chain length
      stressstd=stressstd+(  stress - (data(:,5))    ).^2; % chain length
    end
  end
  ystd=sqrt(ystd/Nseeds);
  Lstd=sqrt(Lstd/Nseeds);
  stressstd=sqrt(stressstd/Nseeds);
  if Nseeds<2
    ysem=zeros(size(ystd));
    Lsem=ysem;
  else 
    ysem=ystd/sqrt(Nseeds-1);
    Lsem=Lstd/sqrt(Nseeds-1);
    stresssem=stressstd/sqrt(Nseeds-1);
  end

  figure
  plot(log10(xdata), log10(ymean)); hold on
  plot(log10(xdata), 0.5*log10(xdata)-1.15)

  ifp=fopen(sprintf('%s/timeprogress_averaged.out', pth), 'w');
  for i=1:length(xdata)
    fprintf(ifp, '%12e %12e %12e %12e %12e %12e %12e %12e\n', xdata(i), ymean(i), ysem(i), Ldata(i), Lsem(i), stress(i), stressstd(i), stresssem(i));
  end
  fclose(ifp);
end
