clear all; 
close all; 
clc
% cd('/home/jschoormans/lood_storage/divi/Temp/Lukas/reconstruction/test-data/spiralshape')
cd('L:\basic\divi\Temp\Lukas\reconstruction\test-data\spiralshape')
file='cs_12012017_1913584_9_2_wipk64r6mode1np50t3flip1senseV4.raw'
% file='cs_12012017_1905246_8_2_wipk64r6mode1np50t3flip0senseV4.raw'
%% copied from lukas' code'
mrecon=MRecon(file)
mrecon.Parameter.Recon.CoilCombination = 'pc';
% cardiac binning
mrecon.Parameter.Cardiac.RetroHoleInterpolation = 'No';
mrecon.Parameter.Cardiac.Synchronization = 'Retrospective';
mrecon.Parameter.Cardiac.RetroPhases = 24;
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
mrecon.K2IM
%%
K=mrecon.Data(32,:,:,:,:,:,:,:,:,:,:,:);
size(K)
% remove stupid checkerboard pattern
che=create_checkerboard([1,64,64]);
K=bsxfun(@times,K,che);
%%

kspace=squeeze(K);
size(kspace)
%%
a = sum(sum(kspace,4),5)./sum(sum(kspace~=0,4),5);
sens=bart('ecalib -r 20 -m1',permute(a,[4 1 2 3]));
%%
close all;
params=params_init();
params.Lg=1;
params.inspectLg=false
params.L3=5;
params.L4=4;
params.sparsity_transform='TV';
params.Imref=[];
params.x=40;
params.y=36;
params.lambda=2e-4
% sens(sens==0)=1e-2;

params.increase_penalty_parameters=false
params.G.precon=false;
params.G.maxiter=10
params.niter=5
P_recon=LRT_recon(kspace,squeeze(sens),params);
%% visualize recon %IS DIFFERENT THAN IN RECON??
figure(1000); immontage4D(squeeze(abs(P_recon)),[]);
figure(1001); immontage4D(squeeze(angle(P_recon)),[-pi pi]);
figure(1002); imshow(angle(P_recon(:,:,1,1,1)),[-pi pi]);

