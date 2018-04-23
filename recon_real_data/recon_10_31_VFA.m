% recon 28 july
clear all; close all; clc

% cd('/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/LRT/Low_Rank_2017_09_21');
cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_10_31')

%%
% MR=MRecon('lr_21092017_1843595_6_2_wipsc23vfat2preptenr18dynsV4.raw')  %OIL
MR=MRecon('lr_31102017_1946079_2_2_wip_brin-vfa-t2prep-lowrankV4.raw') %BRAIN

MR.Parameter.Parameter2Read.typ = 1;
% load data
disp('readdata')
MR.ReadData;
MR.RandomPhaseCorrection;
disp('corrections...')
MR.RemoveOversampling;
MR.PDACorrection; %???
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;

disp('sortdata')

MR.Parameter.Labels.Index.aver=MR.Parameter.Labels.Index.rf;
MR.Parameter.Recon.ImmediateAveraging='No'
MR.SortData;
%%
K=MR.Data{1};
K=ifftshift(ifft(K,[],1),1); size(K)
K=K(80,:,:,:,:,:,:,:,:,:,:,:);
% remove stupid checkerboard pattern
che=create_checkerboard([1,size(K,2),size(K,3)]);
K=bsxfun(@times,K,che);
a = sum(K(:,:,:,:,1,:,1,1,1,1),6)./sum(K(:,:,:,:,1,:,1,1,1,1)~=0,6);
sens=bart('ecalib -S -m1',a);
% sens=ones(size(kspace));
kspace=squeeze(K);

if ndims(kspace)==4 % if one chan
   kspace=permute(kspace,[1 2 5 3 4]); 
end

kspace=kspace(:,:,:,:,:);
a = sum(sum(kspace,4),5)./sum(sum(kspace~=0,4),5);
sens=bart('ecalib -S -m1',permute(a,[4 1 2 3]));
sens=sens+1e-7; % no zero vals in sense maps...

sens=sens;
close all;
params=params_init();
params.Lg=6;
params.inspectLg=false
params.L3=5;
params.L4=5;
params.sparsity_transform='TV';
params.Imref=[];
params.x=17;
params.y=62;
params.lambda=4e-2
params.mu=1e-4;
% sens(sens==0)=1e-2;

params.increase_penalty_parameters=false
params.G.precon=false;
params.G.maxiter=20
params.niter=20;

coil=2
P_recon=LRT_recon(kspace(:,:,:,:,:),squeeze(sens(:,:,:,:)),params);
% visualize recon
figure(1000); immontage4D(squeeze(abs(P_recon)),[]);
%%
x1=30
x2=73
x3=40
y1=32
y2=42
y3=55

figure(200);clf;subplot(221); hold on 
imshow(squeeze(abs(P_recon(:,:,1,5,1))),[]);
plot(y1,x1,'r+')
plot(y2,x2,'b+')
plot(y3,x3,'g+'); hold off
title('first echo, lowest TE');

subplot(223); hold on 
imshow(squeeze(abs(P_recon(:,:,1,1,20))),[]);
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
T2PREP_data=((P_recon(:,:,1,:,1)));
size(T2PREP_data)
T2PREP=[70:-10:30];
[T2_mono_all_allslice, T2_mono_all_rgb_allslice, rsquare_mono_all_allslice] = T2_fitting(T2PREP_data, T2PREP,0.02)

figure(3000); imshow(T2_mono_all_allslice,[10 120]); colormap parula 
colorbar; title('T2 map LRT')
figure(3001); hist(T2_mono_all_allslice((T2_mono_all_allslice>0)),200); title('T2 values')