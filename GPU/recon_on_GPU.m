% test recon 
cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2018_01_31_VFA');
clear
% load('tesdata.mat')
load('testdata_5vc.mat')
params.lambda=1e-1;
params.GPU=1; 
params.subspacedim1=2;
params.G.maxiter=40
params.Lg=6;
params.L3=4;
params.niter=5;
params.alpha=20
P_recon=LRT_recon(kspaceinput,squeeze(sens),params);
%%
figure(1000); immontage4D(squeeze(abs(P_recon)),[]);

%%
 
x1=85
x2=115
x3=155
x4=180; 
y1=18
y2=20
y3=20
y4=60; 

figure(200);clf;subplot(221); hold on 
imshow(squeeze(abs(P_recon(:,:,1,6,1))),[]);
plot(y1,x1,'r+')
plot(y2,x2,'b+')
plot(y3,x3,'g+');
plot(y4,x4,'y+');
hold off
title('first echo, lowest TE');

subplot(223); hold on 
imshow(squeeze(abs(P_recon(:,:,1,1,end))),[]);
plot(y1,x1,'r+')
plot(y2,x2,'b+')
plot(y3,x3,'g+'); 
plot(y4,x4,'y+'); hold off
title('last echo, highest TE');
subplot(222)
hold on
plot([300,200,100,60,40,20,11],squeeze(abs(P_recon(x1,y1,1,:,1))),'r')
plot([300,200,100,60,40,20,11],squeeze(abs(P_recon(x2,y2,1,:,1))),'b')
plot([300,200,100,60,40,20,11],squeeze(abs(P_recon(x3,y3,1,:,1))),'g')
plot([300,200,100,60,40,20,11],squeeze(abs(P_recon(x4,y4,1,:,1))),'y')

hold off
title('T2prep dimension')
ylabel('intensity'); xlabel('TE')

subplot(224)
hold on
plot(squeeze(abs(P_recon(x1,y1,1,6,:))),'r')
plot(squeeze(abs(P_recon(x2,y2,1,6,:))),'b')
plot(squeeze(abs(P_recon(x3,y3,1,6,:))),'g')
plot(squeeze(abs(P_recon(x4,y4,1,6,:))),'y')

hold off
title('echo dimension')
ylabel('intensity'); xlabel('echo number')

%% T2 FITTING 
for i=5
T2PREP_data=(abs(gather(P_recon(:,:,1,:,i))));
size(T2PREP_data)
T2PREP=[300,200,100,70,50,20,10];
[T2_mono_all_allslice, T2_mono_all_rgb_allslice, rsquare_mono_all_allslice] = T2_fitting(T2PREP_data, T2PREP,0.02)
% T2(:,:,i)=T2_mono_all_allslice; 
end
%%
figure(4000); imshow(median(T2_mono_all_allslice,3),[10 400]); colormap parula 
colorbar; title('T2 map LRT')
figure(3001); hist(T2_mono_all_allslice((T2_mono_all_allslice>0)),200); title('T2 values')
 
 
