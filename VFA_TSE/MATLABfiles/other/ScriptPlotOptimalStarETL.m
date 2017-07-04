% Script to investigate optimal ETL for CUBE using the Extended Phase Graph 
% algorithm as described by Busse (2006)

%% Set parameters
T1 = 1200;
T2=35;
TE = 5;
etl = 67;
etlrange=32:80;
% approximation of neck T2 values
ensemble=[20 20 20 30 30 30 30 30 30 30 30 30 30 40 40 40 40 40 40 40 50 50 60 60 60 70 70];
sig=zeros(size(etlrange,2),size(ensemble,2),max(etlrange));



%% do calculations for all different ETLs
for etl=etlrange
    display(['etl = ' num2str(etl) ' of ' num2str(max(etlrange))])
    Star=0.03;
    k=1;
    j=1;
    s=1;
    lastangle=0;
    
%    first optimize Target Sig with given settings 
   while lastangle < pi && s < 30

        Star = Star + 0.025;
        angles = CUBEangles(T1,T2,TE,Star,etl);
        lastangle = angles(end);
        k=k+1;
   end
    
    Star = Star-0.025;
    angles = CUBEangles(T1,T2,TE,Star,etl);
    lastangle = angles(end);

    while lastangle < pi && k < 30

        Star = Star + 0.01;
        angles = CUBEangles(T1,T2,TE,Star,etl);
        lastangle = angles(end);
        k=k+1;
    end

    Star = Star-0.01;
    angles = CUBEangles(T1,T2,TE,Star,etl);
    lastangle = angles(end);

    while lastangle < pi && j < 30
        Star = Star + 0.001;
        angles = CUBEangles(T1,T2,TE,Star,etl);
        lastangle = angles(end);
        j=j+1;
    end

    T2 = T2;
    angles = CUBEangles(T1,T2,TE,Star,etl);
    lastangle = angles(end);
    targg(T2)=Star;
    
% Simulate signal evolution    
    for T2t=1:size(ensemble,2)
        FM=CUBEmagn(T1,ensemble(T2t),TE,angles,etl,1);
        sig(etl,T2t,1:etl)=FM(round(size(FM,1)/2),:);
        sigmean(etl,T2t)=mean(sig(etl,T2t,5:etl));
        sigstd(etl,T2t)=std(sig(etl,T2t,5:etl));
    end
    
end


%% plot results
meansig=mean(sigmean,2);
stdsig=std(sigstd,0,2);

meancov=mean(sigstd./sigmean,2);
stdcov=std(sigstd./sigmean,0,2);

scantime=10176./etlrange;
snrtime=meansig(etlrange)./sqrt(scantime)';

scrsz=get(0,'screensize');
a=figure('Name','Overlay with T1 map (milliseconds)', 'position', [scrsz(3)/9 scrsz(4)/15 scrsz(3)/1.3 scrsz(4)/1.4]);

subplot(1,3,3)
hold on
[ax,H1,H2]=plotyy(etlrange,meansig(etlrange),etlrange,scantime);
set(H1,'LineStyle','--','LineWidth',3)
set(H2,'LineWidth',3)
legend('SNR', 'Scan Time','test')
hold on
errorbar(etlrange,meansig(etlrange),stdsig(etlrange),'b')
ylabel(ax(1),['mean signal of neck spin ensemble'],'fontsize', 8)
ylabel(ax(2),['scan time per CUBE (s)'],'fontsize', 8)
xlabel('Echo train length','fontsize', 8)
title('Signal & Scan Time','fontsize',8)

subplot(1,3,1)
errorbar(etlrange,meancov(etlrange),stdcov(etlrange))
hold on
plot(etlrange,meancov(etlrange),'b','LineWidth',3)
xlabel('echo train length','fontsize', 8)
ylabel(['mean and std of std/mean of signal of neck spin ensemble'],'fontsize', 8)
title('Signal variation (stability)','fontsize',8)

subplot(1,3,2)
plot(etlrange,snrtime,'k','LineWidth', 3)
xlabel('echo train length','fontsize', 8)
ylabel(['snr / sqrt(scan time)'],'fontsize', 8)
axis([5 max(etlrange)+10 0.01 max(snrtime)+0.001])
title('SNR / sqrt(Scan Time)','fontsize',8)

