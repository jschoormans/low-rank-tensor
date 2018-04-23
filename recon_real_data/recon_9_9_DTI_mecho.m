clear; close all; clc
if ispc
    cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_08_13')
    cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_09_05\2017_09_05\lr_23248')
%   addpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations\tensor\low-rank-tensor'))
else
    cd(['/home/',getenv('USER'),'/lood_storage/divi/Ima/parrec/Kerry/LRT_Data/2017_09_09'])
end
%%
clear MR
MR=MRecon('lr_09092017_1919570_3_2_wipdtimechocslowresV4.raw')

DTI=0;
MR.Parameter.Labels.Index.aver=(MR.Parameter.Labels.Index.rf);
MR.Parameter.Parameter2Read.aver=[0:max(MR.Parameter.Labels.Index.aver)].'
MR.Parameter.Recon.ArrayCompression='No';
MR.Parameter.Recon.ACNrVirtualChannels=6;
MR.Parameter.Parameter2Read.typ = 1;
MR.Parameter.Recon.ImmediateAveraging='No';
% MR.Parameter.Parameter2Read.chan=[2;3]
% load data
disp('readdata')
MR.ReadData;
MR.RandomPhaseCorrection;
disp('corrections...')
MR.RemoveOversampling;
MR.PDACorrection; %???
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;
%% MRecon sort
disp('sortdata')
MR.SortData;
K=MR.Data;
Ktemp=squeeze(MR.Data);
K=squeeze(K);
K=fftshift(ifft(ifftshift(K,1),[],1),1);
size(K)

recon_x_loc = 46;
K_1x = K(recon_x_loc,:,:,:,:,:);
Kcc=bart('cc -p5',permute(K_1x,[1 2 3 4 7 8 9 10 5 6]));
Kcc=permute(Kcc,[1 2 3 4 9 10 5 6 7 8]);
size(Kcc)
%% self sort: small memeory allocation
disp('sortdata')
%ifft in readout direction + select slice 
MR.Data=fftshift(ifft(ifftshift(MR.Data,1),[],1),1);

recon_x_loc = 46;
MR.Data=MR.Data(recon_x_loc,:);
K= sortArray(MR);
Kcc=bart('cc -p5',permute(K,[1 2 3 4 7 8 9 10 5 6]));
Kcc=permute(Kcc,[1 2 3 4 9 10 5 6 7 8]);
size(Kcc)

%% all channel recon
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

%% one channel recon 
temp=squeeze(Kcc(1,1:end,1:end,:,:,:,:,:,:,:,:,:));

kspace = temp(:,:,1,:,:);

% remove stupid checkerboard pattern
che=create_checkerboard([size(kspace,1),size(kspace,2)]);
kspace=bsxfun(@times,kspace,che);
sens = ones(size(kspace,1), size(kspace, 2));
%%
params=params_init();
params.L3=5;
params.L4=3;
params.subspacedim1=1;
params.subspacedim2=1; 
params.scaleksp=false; 

params.Lg=2;
params.inspectLg=false;
params.sparsity_transform='TV';
params.Imref=[];
params.x=26;
params.y=34;
params.mu=0.1e2;
params.lambda=5e-3;
% sens(sens==0)=1e-2;

params.niter=5; 
params.increase_penalty_parameters=false;
params.G.precon=true;
params.G.maxiter=10;

P_recon=LRT_recon(kspace,squeeze(sens),params);
%% visualize recon
figure(1000); immontage4D(squeeze(abs(P_recon)),[]);

