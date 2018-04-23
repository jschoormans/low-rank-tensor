clear all;
close all; 
clc
% cd('/home/jschoormans/lood_storage/divi/Temp/Lukas/reconstruction/test-data/spiralshape')
% cd('L:\basic\divi\Temp\Lukas\reconstruction\test-data\spiralshape')
cd('/home/barunderkamp/lood_storage/divi/Temp/Lukas/reconstruction/test-data/spiralshape')
file='cs_12012017_1913584_9_2_wipk64r6mode1np50t3flip1senseV4.raw'
% cd('/home/barunderkamp/lood_storage/divi/Projects/4dflow_lrt/data/2017_10_10_FlowPhantom_Bobby/2017_10_10/Fl_127943')
% file='fl_10102017_1931244_11_2_wip_2dfully_30card_fov_80_nsa4_afprotocolV4.raw'
% cd('/home/barunderkamp/lood_storage/divi/Projects/4dflow_lrt/data/2017_09_27_1stScanday_FlowPhantom_Bobby/2017_09_27/Bo_112719')
% file='bo_27092017_1906118_6_2_wip_3dx4_30cardV4.raw'
% file='cs_12012017_1905246_8_2_wipk64r6mode1np50t3flip0senseV4.raw'
%% copied from lukas' code'
mrecon=MRecon(file)
mrecon.Parameter.Recon.CoilCombination = 'pc';
% cardiac binning
mrecon.Parameter.Cardiac.RetroHoleInterpolation = 'No';
mrecon.Parameter.Cardiac.Synchronization = 'Retrospective';
% mrecon.Parameter.Cardiac.RetroPhases = 24;
mrecon.Parameter.Cardiac.RetroPhases = 30;
mrecon.Parameter.Parameter2Read.typ = 1;
mrecon.Parameter.Parameter2Read.Update;
% load data
mrecon.ReadData;
mrecon.RandomPhaseCorrection;
mrecon.RemoveOversampling;
mrecon.PDACorrection;
mrecon.DcOffsetCorrection;
mrecon.MeasPhaseCorrection;
mrecon.SortData;

% undersamplingmask = makemask(size(mrecon.Data,1),size(mrecon.Data,2));
% % apply retrospective undersampling
% arraymasked = multiplymask(mrecon.Data,undersamplingmask);
% mrecon.Data = arraymasked;

% Version from Jasper
% undersampling_mask=bart(' poisson -Y160 -Z160 -y10 -z10 -C15 -v -e');
% undersampling_mask = permute(undersampling_mask,[2 3 1]);
% mrecon.Data = bsxfun(@times,mrecon.Data,undersampling_mask);

mrecon.K2IM

%%
K=mrecon.Data(32,:,:,:,:,:,:,:,:,:,:,:);
% K=mrecon.Data(:,:,16,:,:,:,:,:,:,:,:,:);

size(K)
% remove stupid checkerboard pattern
% che=create_checkerboard([160,160,1]);
che=create_checkerboard([1,size(K,2),size(K,3)]); % 

K=bsxfun(@times,K,che);
%%

kspace=squeeze(K);
size(kspace)
%%
a = sum(sum(kspace,4),5)./sum(sum(kspace~=0,4),5);
sens=bart('ecalib -r 20 -m1',permute(a,[4 1 2 3]));
%%
sens=sens+1e-4;
close all;
params=params_init();
params.Lg=4;
params.inspectLg=false
params.L3=4;
params.L4=4;
params.sparsity_transform='TV';
params.Imref=[];
params.x=40;
params.y=35;
% params.x=80;
% params.y=16;
% params.y = 80;
params.lambda=8e-1
params.mu=2e0
% sens(sens==0)=1e-2;

params.increase_penalty_parameters=false
params.G.precon=false;
params.G.maxiter=8
params.niter=20
params.subspacedim1=2
params.subspacedim2=8

P_recon=LRT_recon(kspace,squeeze(sens),params);
% visualize recon %IS DIFFERENT THAN IN RECON??
figure(1000); immontage4D(squeeze(abs(P_recon)),[]);
figure(1001); immontage4D(squeeze(angle(P_recon)),[-pi pi]);
figure(1002); imshow(angle(P_recon(:,:,1,1,1)),[-pi pi]);
%
cplkdiff=squeeze((P_recon(:,:,:,:,1))./(P_recon(:,:,:,:,4)));

  figure(1003); immontage4D(reshape(angle(cplkdiff),[64 64 6 5]),[-pi/8 pi/8])
%   figure(1003); immontage4D(reshape(angle(cplkdiff),[160 160 6 5]),[-pi/8 pi/8])

% cplkdiff=squeeze((P_recon(:,:,:,:,1))./(P_recon(:,:,:,:,4)));
%   figure(1003); immontage4D(reshape(angle(cplkdiff),[160 160 6 5]),[-pi/8 pi/8])


%%
% BART RECON 
kspacebart=permute(kspace,[6 1 2 3 7 4 5]);
size(kspacebart)
a = sum(sum(kspacebart,5),6)./sum(sum(kspacebart~=0,5),6);
sensbart=bart('ecalib  -r 20 -m1',a);
size(sensbart)
reconbart=bart('pics -RT:7:0:0.01 -i10 -d5',kspacebart,sens);
figure(2000); immontage4D(abs(squeeze(reconbart)))
figure(2001); immontage4D(angle(squeeze(reconbart)),[-pi,pi])



