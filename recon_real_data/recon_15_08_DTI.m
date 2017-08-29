clear; close all; clc
if ispc
    cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_08_15\')
else
    cd('/home/qzhang/lood_storage/divi/Ima/parrec/Jasper/LRT/Low_Rank_2017_08_15')
    addpath(genpath('/opt/amc/bart/')); vars;
end

MR=MRecon('lr_15082017_2131128_7_2_wip_dti-t2prep_csV4.raw')
%%
DTI=1;

if(~DTI)
    MR.Parameter.Labels.Index.aver=(MR.Parameter.Labels.Index.rf);
    MR.Parameter.Parameter2Read.aver=[0:max(MR.Parameter.Labels.Index.aver)].';
end
MR.Parameter.Recon.ArrayCompression='Yes';
MR.Parameter.Recon.ACNrVirtualChannels=4;
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
MR.SortData;


K=MR.Data;
Ktemp=squeeze(MR.Data); 
clear MR
K=squeeze(K);
if(DTI)
    K = sum(K, 6); %as the nsa label for DTI data set is not correct    
end
K=fftshift(ifft(ifftshift(K,1),[],1),1); 

size(K)

%%
kspace=K(32,1:end,1:end,:,:,:,:,:,:,:,:,:);
imshow(squeeze(abs(kspace(1,:,:,1,1))),[0 1e-2])
size(kspace)

% remove stupid checkerboard pattern
che=create_checkerboard([1,size(kspace,2),size(kspace,3)]);
kspace=bsxfun(@times,kspace,che);
kspace=squeeze(kspace);

if DTI
    par_dim1 = 9; par_dim2 = 6;
    du = zeros(size(kspace,1),size(kspace,2),size(kspace,3),par_dim1,par_dim2);
    du_temp = zeros(size(Ktemp,1),size(Ktemp,2),size(Ktemp,3),size(Ktemp,4),par_dim1,par_dim2);
    for ii = 1:par_dim1
        for jj = 1:par_dim2
            du(:,:,:,ii,jj) = kspace(:,:,:,((jj-1) * par_dim1 +ii));
            du_temp(:,:,:,:,ii,jj) = Ktemp(:,:,:,:,((jj-1) * par_dim1 +ii));
        end
    end
    kspace=du;
    Ktemp = du_temp;
end
figure(111); immontage4D(squeeze(abs(kspace(:,:,1,:,:))>0),[]);
%%
% a = sum(sum(kspace(:,:,:,:,:),4),5)./sum(sum(kspace(:,:,:,:,:)~=0,4),5);
a=kspace(:,:,:,1,6);  %highest SNR volume
sens=bart('ecalib -S -m1 -c0.1',permute(a,[4 1 2 3]));
% sens=sens+1e-7; % no zero vals in sense maps...

figure(112)
subplot(311)
immontage4D(abs(sens),[])
subplot(312)
immontage4D(real(sens),[])
subplot(313)
immontage4D(angle(sens),[-pi pi])

%% 3,4, ok, 2 too noisy
coilnr=[3]
sens_onecoil=sens(:,:,:,coilnr); 
kspace_onecoil=kspace(:,:,coilnr,:,:);

params=params_init();
params.L3=3;
params.L4=4;
params.subspacedim1=1;
params.subspacedim2=6; 
[nav_estimate_1,nav_estimate_2,eigenvals_1,eigenvals_2]= subspace_estimate_3D(Ktemp(:,:,:,2,:,:),params);

params.nav_estimate_1=nav_estimate_1;
params.nav_estimate_2=nav_estimate_2;
params.eigenvals_1=eigenvals_1;
params.eigenvals_2=eigenvals_2;

params.Lg=2;
params.inspectLg=true;
params.sparsity_transform='TV';
params.Imref=[];
params.x=40;
params.y=65;

params.mu=1e1;
params.lambda=5e-3;
params.automu=1; 
params.autolambda=1;

params.scaleksp=0
params.niter=5; 
params.increase_penalty_parameters=false;
params.G.precon=false;
params.G.maxiter=5;
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
