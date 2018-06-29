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

%%
mrecon.K2IM
%%
K=mrecon.Data(32,:,:,:,:,:,:,:,:,:,:,:);
size(K)
% remove stupid checkerboard pattern
che=create_checkerboard([1,64,64]);
K=bsxfun(@times,K,che);

a = sum(K(:,:,:,:,1,:,1,1,1,1),6)./sum(K(:,:,:,:,1,:,1,1,1,1)~=0,6);
sens=bart('ecalib -r 20 -m1',a);

%%
%%
params=params_init();
params.L3=4;
params.L4=4;
params.subspacedim1=2
params.subspacedim2=2

params.scaleksp=true; 

params.Lg=12;
params.inspectLg=false;
params.sparsity_transform='TV';
params.Imref=[];
params.x=17;
params.y=39;
params.mu=4;
params.lambda=0.5;
% sens(sens==0)=1e-2;

params.niter=3; 
params.increase_penalty_parameters=false;
params.G.maxiter=50;
params.C.maxiter=50;
params.alpha=10
params.mu=10;


params.TVoption=2

P_recon=LRT_recon_test(squeeze(K),squeeze(sens),params);

imagine(squeeze(angle(P_recon)))