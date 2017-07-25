clear all; close all; clc
cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_07_24\2017_07_24\lr_2407')
addpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations\tensor\low-rank-tensor'))
MR=MRecon('lr_24072017_1920025_3_2_wipvfat2preplowrankV4.raw')

MR.Parameter.Labels.Index.aver=(MR.Parameter.Labels.Index.rf);
MR.Parameter.Labels.NumberOfEchoes=60
MR.Parameter.Recon.RemoveMOversampling='No'
MR.Parameter.Recon.RemovePOversampling='No'
MR.Parameter.Recon.ArrayCompression='Yes'
MR.Parameter.Recon.ACNrVirtualChannels=4
MR.Parameter.Parameter2Read.typ = 1;
MR.Parameter.Recon.ImmediateAveraging='No'
MR.Parameter.Parameter2Read.Update;

% load data
MR.ReadData;
MR.RandomPhaseCorrection;
MR.RemoveOversampling;
MR.PDACorrection;
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;
MR.SortData;
MR.K2IM;

K=MR.Data(70,:,:,:,:,:,:,:,:,:,:,:);
che=create_checkerboard([1,size(K,2),size(K,3)]);
K=bsxfun(@times,K,che);
a = sum(K(:,:,:,:,1,:,1,1,1,1),6)./sum(K(:,:,:,:,1,:,1,1,1,1)~=0,6);
sens=bart('ecalib -S -m1',a);
%%
params=params_init();
params.Lg=2;
params.L3=4; 
params.L4=6;
params.mu=1e4
params.sparsity_transform='TV'
LRT_recon(squeeze(K),squeeze(sens),params)