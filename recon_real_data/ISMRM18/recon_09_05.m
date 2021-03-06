clear; close all; clc
if ispc
      cd('L:\basic\divi\Ima\parrec\Kerry\LRT_Data\2017_09_05_knee')

%     addpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations\tensor\low-rank-tensor'))
else
    cd(['/home/',getenv('USER'),'/lood_storage/divi/Ima/parrec/Jasper/LRT/Low_Rank_2017_09_05'])
    addpath(genpath('/opt/amc/bart/')); vars;
end
%%
clear MR
MR=MRecon('lr_05092017_2032554_12_2_wipvfat2prepcstransverseV4.raw')
DTI=0;

MR.Parameter.Labels.Index.aver=(MR.Parameter.Labels.Index.rf);
MR.Parameter.Parameter2Read.aver=[0:max(MR.Parameter.Labels.Index.aver)].'
% MR.Parameter.Parameter2Read.aver=[0:6].'

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
MR.Data=MR.Data(93,:); %maybe 93??
K= sortArray(MR);
Kcc=bart('cc -p4',permute(K,[1 2 3 4 7 8 9 10 5 6]));
Kcc=permute(Kcc,[1 2 3 4 9 10 5 6 7 8]);
size(Kcc)
%%
kspace=Kcc(1,1:end,1:end,:,:,:,:,:,:,:,:,:);
% kspace=K(80,2:end,2:end,:,:,:,:,:,:,:,:,:);
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
params.subspacedim1=1;
params.subspacedim2=5; 
params.scaleksp=false; 

params.Lg=3;
params.inspectLg=false;
params.sparsity_transform='TV';
params.Imref=[];
params.x=20;
params.y=45;
params.mu=0.1e2;
params.automu=1;
params.autolambda=1;
params.lambda=5e-3;
% sens(sens==0)=1e-2;

params.niter=2; 
params.increase_penalty_parameters=false;
params.G.precon=true;
params.G.maxiter=10;

P_recon=LRT_recon(kspace,squeeze(sens),params);
%% visualize recon
figure(1000); immontage4D(squeeze(abs(P_recon)),[]);

