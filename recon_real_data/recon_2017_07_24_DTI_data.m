clear all; close all; clc
%% LOAD DATA
if ispc
%     cd('L:\basic\divi\Ima\parrec\Jasper\Low_Rank_2017_07_11')
    cd('L:\basic\divi\Ima\parrec\Jasper\Low_Rank_2017_07_17')
    addpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations\tensor\low-rank-tensor'))
else
    cd(['/home/',getenv('USER'),'/lood_storage/divi/Ima/parrec/Kerry/Data/2017_07_24_LRT_DTI_T2prep'])  
end

% MR = MRecon('lo_12072017_1926046_17_2_wipvfat2preplowrankV4.raw');
% MR=MRecon('lo_17072017_1603264_2_2_wipvfat2preplowrankV4.raw')
% MR=MRecon('lo_17072017_1635333_4_2_wipvfat2preplowrankV4.raw')
MR=MRecon('lr_24072017_2013158_9_2_wipvfat2preplowrankV4.raw')

%     

MR.Parameter.Recon.RemoveMOversampling='No';
MR.Parameter.Recon.RemovePOversampling='No';
MR.Parameter.Recon.ArrayCompression='Yes';
MR.Parameter.Recon.ACNrVirtualChannels=4;

MR.Parameter.Parameter2Read.typ = 1;
MR.Parameter.Recon.ImmediateAveraging='No'

% MR.Parameter.Recon.CoilCombination = 'pc';
% MR.Parameter.Cardiac.RetroHoleInterpolation = 'No';
% MR_data.Parameter.Cardiac.Synchronization = 'Yes';
% MR.Parameter.Cardiac.RetroPhases = 24;

MR.Parameter.Parameter2Read.Update;

% load data
disp('readdata')
MR.ReadData;
MR.RandomPhaseCorrection;
disp('corrections...')
MR.RemoveOversampling;
MR.PDACorrection;
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;

disp('sortdata')
MR.SortData;
MR.K2IM;
%% Preprocess k space
% K=bart('cc -p3',MR.Data);
data = double(MR.Data);
K=MR.Data(124,:,:,:,:,:,:,:,:,:,:,:);
K = sum(K, 12);
% remove stupid checkerboard pattern
che=create_checkerboard([1,size(K,2),size(K,3)]);
K=bsxfun(@times,K,che);

% a = sum(K(:,:,:,:,1,:,1,1,1,1),6)./sum(K(:,:,:,:,1,:,1,1,1,1)~=0,6);
a = K(:,:,:,:,1);
sens=bart('ecalib -S -m1',a);
du_temp=squeeze(K);

par_dim1 = 9; par_dim2 = 6;
du = zeros(size(du_temp,1),size(du_temp,2),size(du_temp,3),par_dim1,par_dim2);
for ii = 1:par_dim1
    for jj = 1:par_dim2
        du(:,:,:,ii,jj) = du_temp(:,:,:,((jj-1) * par_dim1 +ii));
    end
end

%% RECON

params = params_init;
params.L3 = 2; 
params.L4 = 3;
params.Lg = 2
LRT_recon(du,squeeze(sens),params);

