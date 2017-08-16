clear; close all; clc
if ispc
    cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_08_15\')
else
    cd('/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/LRT/Low_Rank_2017_08_13')
    addpath(genpath('/opt/amc/bart/')); vars;
end

MR=MRecon('lr_15082017_2116431_6_2_wip_vfa-t2prep_csV4.raw')
%%
DTI=0;

MR.Parameter.Labels.Index.aver=(MR.Parameter.Labels.Index.rf);
MR.Parameter.Parameter2Read.aver=[0:max(MR.Parameter.Labels.Index.aver)].'
MR.Parameter.Recon.ArrayCompression='Yes';
MR.Parameter.Recon.ACNrVirtualChannels=4;
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
K=fftshift(ifft(ifftshift(K,1),[],1),1); 
size(K)
Ktemp=squeeze(MR.Data); 

%%
kspace=K(60,1:end,1:end,:,:,:,:,:,:,:,:,:);
imshow(squeeze(abs(kspace(1,:,:,1,1))),[0 1e-2])
size(kspace)


% remove stupid checkerboard pattern
che=create_checkerboard([1,size(kspace,2),size(kspace,3)]);
kspace=bsxfun(@times,kspace,che);
kspace=squeeze(kspace);

if DTI
par_dim1 = 9; par_dim2 = 6;
du = zeros(size(kspace,1),size(kspace,2),size(kspace,3),par_dim1,par_dim2);
for ii = 1:par_dim1
    for jj = 1:par_dim2
        du(:,:,:,ii,jj) = kspace(:,:,:,((jj-1) * par_dim1 +ii));
    end
end
kspace=du; 
end
%%
% a = sum(sum(kspace(:,:,:,:,:),4),5)./sum(sum(kspace(:,:,:,:,:)~=0,4),5);
a=kspace(:,:,:,5,1);
sens=bart('ecalib -S -r25 -m1 -c0.015',permute(a,[4 1 2 3]));
% sens=sens+1e-7; % no zero vals in sense maps...

figure(112)
subplot(311)
immontage4D(abs(sens),[])
subplot(312)
immontage4D(real(sens),[])
subplot(313)
immontage4D(angle(sens),[-pi pi])

%%

sens_onecoil=sens(:,:,:,[3]); 
kspace_onecoil=kspace(:,:,[3],:,:);

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

params.Lg=1;
params.inspectLg=false;
params.sparsity_transform='TV';
params.Imref=[];
params.x=50;
params.y=55;

params.mu=1e2;
params.lambda=10e0-3;
params.automu=1; 
params.autolambda=1; 

% sens(sens==0)=1e-2;

params.niter=5; 
params.increase_penalty_parameters=false;
params.G.precon=true;
params.G.maxiter=10;

P_recon=LRT_recon(kspace_onecoil,squeeze(sens_onecoil),params);
%% visualize recon
figure(1000); immontage4D(squeeze(abs(P_recon)),[]);
figure(1001);imshow(abs(P_recon(:,:,5,1)).',[])
