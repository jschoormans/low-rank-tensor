% Script to investigate optimal T2 for flip angle calculation for CUBE using 
% the Extended Phase Graph algorithm as described by Busse (2006)

%% Set parameters

T1 = 1200;
TE = 5;
etl = 67;
t2range=10:80;
% approximation of neck T2 values
ensemble=[20 20 20 30 30 30 30 30 30 30 30 30 30 40 40 40 40 40 40 40 50 50 60 60 60 70 70];

clear sig
% sig=zeros(size(t2range,2),size(ensemble,2),max(t2range));

%% do calculations
for T2=t2range
    display(['T2 = ' num2str(T2) ' of ' num2str(max(t2range))])
    Star=0.05;
    k=1;
    kk=1;
    j=1;
    s=1;
    lastangle=0;
    
    %first optimize target signal
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
    
    
    while lastangle < pi && kk < 30
        Star = Star + 0.0001;
        angles = CUBEangles(T1,T2,TE,Star,etl);
        lastangle = angles(end);
        kk=kk+1;
    end

%     T2 = T2;
    Star = Star-0.0001;
    angles = CUBEangles(T1,T2,TE,Star,etl);
    lastangle = angles(end);
    targg(T2)=Star;
    

    
    for T2t=1:size(ensemble,2)
        FM=CUBEmagn(T1,ensemble(T2t),TE,angles,etl,1);
        sig(T2,T2t,:)=FM(round(size(FM,1)/2),:);
    end
    
end

%% plot results
scrsz=get(0,'screensize');
a=figure('Name','Overlay with T1 map (milliseconds)', 'position', [scrsz(3)/9 scrsz(4)/15 scrsz(3)/2 scrsz(4)/3]);
subplot(1,3,1)
% sigdec=sig(:,:,5)./sig(:,:,64);
sigmean=mean(sig(:,:,5:end),3);
sigstd=std(sig(:,:,5:end),0,3);
sigcov=sigstd./sigmean;
meansigcov=mean(sigcov,2);
stdsigcov=std(sigcov,0,2);
% sigdecstd=std(sigdec');
% sigdecmean=mean(sigdec,2);
% sigdecstdmean=sigdecstd./sigdecmean';
%plot sig difference between echo 15 and 60 with cov (std/mean)
errorbar(t2range,meansigcov(t2range),stdsigcov(t2range),'r')
hold on
plot(t2range,meansigcov(t2range),'r','LineWidth', 3)
set(gca,'fontsize',8)
title('Signal variation (stability)','fontsize',8)
axis([5 85 0.05 0.4])
xlabel('T2 used for flip angle calculation (ms)','fontsize', 8)
ylabel(['mean and std of std/mean of signal of neck spin ensemble'],'fontsize', 8)

sigmag=mean(mean(sig(:,:,5:end),3),2);
sigmagstd=std(mean(sig(:,:,5:end),3)');
sigmagstdmean=sigmagstd./sigmag';
subplot(1,3,2)
%plot mean signal level with cov (std/mean)
errorbar(t2range,sigmag(t2range),sigmagstd(t2range),'b')
hold on
plot(t2range,sigmag(t2range),'--b','LineWidth',3)
set(gca,'fontsize',8)
title('SNR','fontsize',8)
axis([5 85 0.1 0.2])
xlabel('T2 used for flip angle calculation (ms)','fontsize', 8)
ylabel(['mean and std of signal of neck spin ensemble'],'fontsize', 8)

subplot(1,3,3)
[ax,H1,H2]=plotyy(t2range,sigmag(t2range),t2range,meansigcov(t2range))
set(H1,'LineStyle','--','LineWidth',3)
set(H2,'LineWidth',3)
ylabel(ax(1),['mean signal of neck spin ensemble'],'fontsize', 8)
ylabel(ax(2),['mean std/mean of signal of neck spin ensemble'],'fontsize', 8)
xlabel('T2 used for flip angle calculation (ms)','fontsize', 8)
legend('SNR', 'Signal variation (stability)')
title('Signal variation (stability) & SNR','fontsize',8)
axis(ax(1),[5 85 0.11 0.18])
axis(ax(2),[5 85 0.05 0.35])




% hold on
% errorbar(10:50,sigdecmean,sigdecstd)
% [ax,h1,h2]=plotyy(10:50,sigdecmean,10:50,sigmag)
% 
% set(gca,'fontsize',8)
% set(ax(1),'fontsize',8)
% set(ax(2),'fontsize',8)
% name = ['SimuCUBE.tif'];
% % box off
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 4])
% print(a, '-r600', '-dtiff', name);









