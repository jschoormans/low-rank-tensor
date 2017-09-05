clear; close all; clc
if ispc
    cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_08_15\')
else
    cd('/home/qzhang/lood_storage/divi/Ima/parrec/Jasper/LRT/Low_Rank_2017_08_15')
end

MR=MRecon('lr_15082017_2116431_6_2_wip_vfa-t2prep_csV4.raw')
%%
DTI=0;

MR.Parameter.Labels.Index.aver=(MR.Parameter.Labels.Index.rf);
MR.Parameter.Parameter2Read.aver=[0:max(MR.Parameter.Labels.Index.aver)].'
MR.Parameter.Recon.ArrayCompression='Yes';
MR.Parameter.Recon.ACNrVirtualChannels=;
MR.Parameter.Parameter2Read.typ = 1;
MR.Parameter.Recon.ImmediateAveraging='No';

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
MR.SortData;


K=MR.Data;
Ktemp=squeeze(MR.Data); 
clear MR
K=squeeze(K);
for i=1:size(K,6); 
K(:,:,:,:,:,i)=fftshift(ifft(ifftshift(squeeze(K(:,:,:,:,:,i)),1),[],1),1); 
end
size(K)

%%
kspace=K(113,1:end,1:end,:,:,:,:,:,:,:,:,:);
imshow(squeeze(abs(kspace(1,:,:,1,1))),[0 1e-2])
size(kspace)

% remove stupid checkerboard pattern
che=create_checkerboard([1,size(kspace,2),size(kspace,3)]);
kspace=bsxfun(@times,kspace,che);
kspace=squeeze(kspace);

%
% a = sum(sum(kspace(:,:,:,:,:),4),5)./sum(sum(kspace(:,:,:,:,:)~=0,4),5);
a=kspace(:,:,:,5,1);
sens=bart('ecalib -S -m1 -c0.1',permute(a,[4 1 2 3]));
% sens=sens+1e-7; % no zero vals in sense maps...

figure(112)
subplot(311)
immontage4D(abs(sens),[])
subplot(312)
immontage4D(real(sens),[])
subplot(313)
immontage4D(angle(sens),[-pi pi])

% quick recon of all 
a2 = sum(sum(kspace(:,:,:,:,:),4),5)./sum(sum(kspace(:,:,:,:,:)~=0,4),5);
a2(isnan(a2))=0;
% im2=ifft(ifft(a2,[],1),[],2);
im2=ifft2c(a2);
im2=sqrt((sum(im2.^2,3)));
im2(isnan(im2))=0; 
figure(113); imshow(abs(im2),[])
%% 3,4, ok, 2 too noisy
coilnr=[3]
sens_onecoil=sens(:,:,:,coilnr); 
kspace_onecoil=kspace(:,:,coilnr,:,:);

params=params_init();
params.L3=3;
params.L4=4;
params.subspacedim1=1;
params.subspacedim2=5; 
[nav_estimate_1,nav_estimate_2,eigenvals_1,eigenvals_2]= subspace_estimate_3D(Ktemp(:,:,:,2,:,:),params);

params.nav_estimate_1=nav_estimate_1;
params.nav_estimate_2=nav_estimate_2;
params.eigenvals_1=eigenvals_1;
params.eigenvals_2=eigenvals_2;

params.Lg=2;
params.inspectLg=false;
params.sparsity_transform='TV';
params.Imref=[];
params.x=40;
params.y=65;

params.mu=1e2;
params.lambda=5e-3;
params.automu=1; 
params.autolambda=1;

params.scaleksp=0
params.niter=5; 
params.increase_penalty_parameters=false;
params.G.precon=false;
params.G.maxiter=15;
params.normalize_sense=1
P_recon=LRT_recon(kspace_onecoil,squeeze((sens_onecoil)),params);
%% visualize recon
figure(1000); immontage4D(squeeze(abs(P_recon)),[]);
figure(1001);imshow(abs(P_recon(:,:,1,1,60)).',[])
figure(1002)
imshow(cat(1,cat(2,...
    abs(P_recon(:,:,1,1,60).'),...
    abs(P_recon(:,:,1,3,60).'),...
    abs(P_recon(:,:,1,5,60).')),...
    cat(2,abs(P_recon(:,:,1,1,60).'),...
    abs(P_recon(:,:,1,3,60).'),...
    abs(P_recon(:,:,1,5,60).'))),[]);

figure(1003);imshow(abs(sqrt(sum(P_recon(:,:,1,3,:).^2,5))).',[])
