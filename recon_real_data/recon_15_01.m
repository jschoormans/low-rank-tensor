clear; close all; clc
cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2018_01_15_VFA')
% recon T2-VFA copied from: recon_10_31_modfrom_5_9.m


%%
clear MR
% MR=MRecon('lr_15012018_1852004_3_2_wip_brin-vfa-t2prep-lowrankV4.raw')
% MR=MRecon('lr_15012018_1930556_3_2_wip_brin-vfa-t2prep-lowrankV4.raw')
% MR=MRecon('lr_15012018_1941075_4_2_wip_brin-vfa-t2prep-lowrankV4.raw')
% MR=MRecon('lr_15012018_1958067_5_2_wip_brin-vfa-t2prep-lowrankV4.raw');
MR=MRecon('lr_15012018_2009293_6_2_wip_brin-vfa-t2prep-lowrankV4.raw')

DTI=0;
MR.Parameter.Labels.Index.aver=(MR.Parameter.Labels.Index.rf);
MR.Parameter.Parameter2Read.aver=[0:max(MR.Parameter.Labels.Index.aver)].'
MR.Parameter.Recon.ArrayCompression='No';
MR.Parameter.Recon.ACNrVirtualChannels=6;
MR.Parameter.Parameter2Read.typ = 1;
MR.Parameter.Recon.ImmediateAveraging='No';
% MR.Parameter.Parameter2Read.chan=[34;35;36;37;44;45]
% load data
disp('readdata')
MR.ReadData;
MR.RandomPhaseCorrection;
disp('corrections...')
MR.RemoveOversampling;
MR.PDACorrection; %???
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;

%%
disp('sortdata')
% MR.SortData;

%ifft in readout direction + select slice 
MR.Data=fftshift(ifft(ifftshift(MR.Data,1),[],1),1);
MR.Data=MR.Data(125,:);
K= sortArray(MR);
Kcc=bart('cc -p5',permute(K,[1 2 3 4 7 8 9 10 5 6]));
Kcc=permute(Kcc,[1 2 3 4 9 10 5 6 7 8]);
size(Kcc)
%%
kspace=Kcc(1,1:end,1:end,:,:,:,:,:,:,:,:,:);
imshow(squeeze(abs(kspace(1,:,:,1,1))),[0 1e-2])
size(kspace)

% remove stupid checkerboard pattern
che=create_checkerboard([1,size(kspace,2),size(kspace,3)]);
kspace=bsxfun(@times,kspace,che);
kspace=squeeze(kspace);

a = sum(sum(kspace(:,:,:,:,:),4),5)./(sum(sum(kspace(:,:,:,:,:)~=0,4),5)+eps);
a2 = (sum(kspace(:,:,:,:,:),5))./((sum(kspace(:,:,:,:,:)~=0,5))+eps);

sens=bart('ecalib -S -r15 -m1',permute(a,[4 1 2 3]));
% sens=sens+1e-7; % no zero vals in sense maps...

figure(112)
subplot(311)
immontage4D(abs(sens),[])
subplot(312)
immontage4D(real(sens),[])
subplot(313)
immontage4D(angle(sens),[-pi pi])

%%
params=params_init();
params.L3=3;
params.L4=3;
params.subspacedim1=50;
params.subspacedim2=5; 
params.scaleksp=false; 

params.Lg=5;
params.inspectLg=false;
params.sparsity_transform='TV';
params.Imref=[];
params.x=45;
params.y=120;
params.mu=0.3e3;
params.lambda=15e-2;
% sens(sens==0)=1e-2;

params.niter=3; 
params.increase_penalty_parameters=false;
params.G.precon=true;
params.G.maxiter=10;

P_recon=LRT_recon(kspace,squeeze(sens),params);
%% visualize recon
figure(1000); immontage4D(squeeze(abs(P_recon)),[]);
 %% TO DO 

 
 
x1=130
x2=173
x3=140
y1=32
y2=33
y3=40

figure(200);clf;subplot(221); hold on 
imshow(squeeze(abs(P_recon(:,:,1,5,1))),[]);
plot(y1,x1,'r+')
plot(y2,x2,'b+')
plot(y3,x3,'g+'); hold off
title('first echo, lowest TE');

subplot(223); hold on 
imshow(squeeze(abs(P_recon(:,:,1,1,50))),[]);
plot(y1,x1,'r+')
plot(y2,x2,'b+')
plot(y3,x3,'g+'); hold off
title('last echo, highest TE');
subplot(222)
hold on
plot([60,50,40,30,20],squeeze(abs(P_recon(x1,y1,1,:,1))),'r')
plot([60,50,40,30,20],squeeze(abs(P_recon(x2,y2,1,:,1))),'b')
plot([60,50,40,30,20],squeeze(abs(P_recon(x3,y3,1,:,1))),'g')
hold off
title('T2prep dimension')
ylabel('intensity'); xlabel('TE')

subplot(224)
hold on
plot(squeeze(abs(P_recon(x1,y1,1,1,:))),'r')
plot(squeeze(abs(P_recon(x2,y2,1,1,:))),'b')
plot(squeeze(abs(P_recon(x3,y3,1,1,:))),'g')
hold off
title('echo dimension')
ylabel('intensity'); xlabel('echo number')

%% T2 FITTING 
for i=1:50
T2PREP_data=((P_recon(:,:,1,:,i)));
size(T2PREP_data)
T2PREP=[70,50,30,20,10];
[T2_mono_all_allslice, T2_mono_all_rgb_allslice, rsquare_mono_all_allslice] = T2_fitting(T2PREP_data, T2PREP,0.02)
T2(:,:,i)=T2_mono_all_allslice; 
end

figure(3000); imshow(mean(T2,3),[10 120]); colormap parula 
colorbar; title('T2 map LRT')
figure(3001); hist(T2_mono_all_allslice((T2_mono_all_allslice>0)),200); title('T2 values')
 
 
 
 
 
 
 
 
 
 
 
 
