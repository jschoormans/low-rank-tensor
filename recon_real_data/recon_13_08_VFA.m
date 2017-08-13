clear; close all; clc
if ispc
    cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_08_13')
%     addpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations\tensor\low-rank-tensor'))
else
    cd('/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/LRT/Low_Rank_2017_08_13')
    addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/tensor/low-rank-tensor'))
%     addpath(genpath('/home/qzhang/lood_storage/divi/Projects/cosart/CS_simulations/tensor/low-rank-tensor'))
    addpath(genpath('/opt/amc/bart/')); vars;
end
MR=MRecon('lr_13082017_1702510_29_2_wip_sc18-vfa-dtiV4.raw')
% MR=MRecon('lr_13082017_1652290_28_2_wip_sc23-vfa-t2prep_iV4.raw')
% MR=MRecon('lr_13082017_1431351_12_2_wip_sc23-vfa-t2prep_iV4.raw')
% MR=MRecon('lr_13082017_1157339_3_2_wip_sc23-vfa-t2prep_tenr18dynsV4.raw')
% MR=MRecon('lr_13082017_1545170_21_2_wip_sc23-vfa-t2prep_iV4.raw')
% MR=MRecon('lr_13082017_1641548_26_2_wip_surveyV4.raw')
%%
DTI=1;

if ~DTI
    MR.Parameter.Labels.Index.aver=(MR.Parameter.Labels.Index.rf);
    MR.Parameter.Parameter2Read.aver=[0:max(MR.Parameter.Labels.Index.aver)].'
    % MR.Parameter.Parameter2Read.aver=[0:60].';
    % MR.Parameter.Labels.NumberOfEchoes=60
    % MR.Parameter.Recon.RemoveMOversampling='No'
    % MR.Parameter.Recon.RemovePOversampling='No'
    MR.Parameter.Recon.ArrayCompression='No';
    MR.Parameter.Recon.ACNrVirtualChannels=6;
    MR.Parameter.Parameter2Read.typ = 1;
    MR.Parameter.Recon.ImmediateAveraging='No';
    MR.Parameter.Parameter2Read.chan=[2,3,4].';
%     MR.Parameter.Parameter2Read.ky=[-50:50].'

else
    %MR.Parameter.Labels.Index.aver=zeros(size(MR.Parameter.Labels.Index.aver));
    %MR.Parameter.Recon.ImmediateAveraging='Yes'
    MR.Parameter.Recon.ArrayCompression='No';
    MR.Parameter.Recon.ACNrVirtualChannels=6;
    MR.Parameter.Parameter2Read.typ = 1;
end

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
%%
%{
params=params_init();
params.L3=4;
params.L4=4;
params.subspacedim1=1;
params.subspacedim2=5; 

[nav_estimate_1,eigenvals_1,nav_estimate_2,eigenvals_2] = full_subspace_estimate_3D(K,params)
%}

%%
K=MR.Data;

K=ifftshift(ifft(K,[],1),1); 
size(K)
K=K(32,:,:,:,:,:,:,:,:,:,:,:);

% remove stupid checkerboard pattern
che=create_checkerboard([1,size(K,2),size(K,3)]);
K=bsxfun(@times,K,che);

kspace=squeeze(K);

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
a = sum(sum(kspace(:,:,:,1,:),4),5)./sum(sum(kspace(:,:,:,1,:)~=0,4),5);
sens=bart('ecalib -S -r20 -m1',permute(a,[4 1 2 3]));
sens=sens+1e-7; % no zero vals in sense maps...
figure(112)
immontage4D(abs(sens),[])
%%
sens=sens;
close all;
params=params_init();
params.Lg=1;
params.inspectLg=false
params.sparsity_transform='TV';
params.Imref=[];
params.x=50;
params.y=25;
params.mu=1e3;
params.lambda=0.5e-3
% sens(sens==0)=1e-2;

params.increase_penalty_parameters=false
params.G.precon=false;
params.G.maxiter=10
P_recon=LRT_recon(kspace,squeeze(sens),params);
%% visualize recon
figure(1000); immontage4D(squeeze(abs(P_recon)),[]);
