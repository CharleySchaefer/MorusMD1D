function getStickerDistribution()
clc; close all;
  fname='seed1/timeprogress.out';
 delim=' ';
  Nheader=1;

  % Get columns with sticker positions
  fid = fopen(fname, 'rt'); 
  tline = fgets(fid);
  fclose(fid);
  header_arr=strsplit(tline, delim);
 ycol=[];
  for i=1:length(header_arr)
  header=header_arr(i){1};
    %(header, (1)
    if( strncmp(header, 'R01',3 ) )
      ycol=[ycol, i];
    end
  end
  ZS=length(ycol);

  % Get all data
  for seed=1:5
seed
  fname=sprintf('seed%d/timeprogress.out', seed);
  data=importdata(fname,delim,Nheader);
  data=data.data;
  
  t=data(:,1);
  R1=data(:,2);
  R2=data(:,3);
if seed ==1  % initialise
  Rmin=min(R1);
  Rmax=max(R2);
  L=Rmax-Rmin;
  Rmin=Rmin-0.5*L;
  Rmax=Rmax+0.5*L;
  NR=200;
  Rspace=linspace(Rmin, Rmax,NR);
  dR=Rspace(2)-Rspace(1);
  Cend=zeros(1,NR);
  CS=zeros(1,NR);
  Rdist=zeros(1,NR);
end

  tstart=2.85e7; counter=0;
  for i=1:length(t)
    if t(i)>tstart
      ind=round( (R1(i)-Rmin)/dR );
      Cend(ind)=Cend(ind)+1;
      ind=round( (R2(i)-Rmin)/dR );
      Cend(ind)=Cend(ind)+1;

      for j=1:ZS
        Rsticker=data(i,ycol(j));
        if j>1
        ind=round( (Rsticker-Rsticker0-Rmin)/dR );
        Rdist(ind)=Rdist(ind)+1; Rsticker0=Rsticker;
          
        end
        ind=round( (Rsticker-Rmin)/dR );
        CS(ind)=CS(ind)+1; Rsticker0=Rsticker;
      end

      counter=counter+1;
    end
  end
  end
  figure
  subplot(2,2,1:2)
  plot(t,R1); hold on
  plot(t,R2); hold on
  subplot(2,2,3)
  plot(Rspace, Cend); hold on
  plot(Rspace, CS); hold on
  subplot(2,2,4)
  plot(Rspace, Rdist); hold on
end
