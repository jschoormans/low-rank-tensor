clear all; close all; clc
if ispc
    cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_07_17')
    addpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations\tensor\low-rank-tensor'))
else
    cd('/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/LRT/Low_Rank_2017_08_07')
    addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/tensor/low-rank-tensor'))
%     cd('/home/qzhang/lood_storage/divi/Ima/parrec/Jasper/Low_Rank_2017_07_24/2017_07_24/lr_2407/')
    addpath(genpath('/home/zhang/lood_storage/divi/Projects/cosart/CS_simulations/tensor/low-rank-tensor'))
    addpath(genpath('/opt/amc/bart/')); vars;
end

MR=MRecon('lr_07082017_1745143_7_2_wip_vfa-t2prep_csV4.raw')

%%
DTI=false;

if ~DTI
    MR.Parameter.Labels.Index.aver=(MR.Parameter.Labels.Index.rf);
    MR.Parameter.Parameter2Read.aver=[0:max(MR.Parameter.Labels.Index.aver)].'
    % MR.Parameter.Parameter2Read.aver=[0:60].';
    % MR.Parameter.Labels.NumberOfEchoes=60
    % MR.Parameter.Recon.RemoveMOversampling='No'
    % MR.Parameter.Recon.RemovePOversampling='No'
    MR.Parameter.Recon.ArrayCompression='No';
    MR.Parameter.Recon.ACNrVirtualChannels=4;
    MR.Parameter.Parameter2Read.typ = 1;
    MR.Parameter.Recon.ImmediateAveraging='No';
%     MR.Parameter.Parameter2Read.chan=[2,3,4].';
%     MR.Parameter.Parameter2Read.ky=[-50:50].'

else
    %MR.Parameter.Labels.Index.aver=zeros(size(MR.Parameter.Labels.Index.aver));
    %MR.Parameter.Recon.ImmediateAveraging='Yes'
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
K=MR.Data;
res=50;
K=K(154-res:155+res,154-res:155+res,:,:,:,:,:,:,:,:,:,:);
size(K)

K=ifftshift(ifft(K,[],1),1); 
disp('LIMITED RESOLUTION')
size(K)
K=K(57,:,:,:,:,:,:,:,:,:,:,:);
K=flip(K,5);
% remove stupid checkerboard pattern
che=create_checkerboard([1,size(K,2),size(K,3)]);
K=bsxfun(@times,K,che);
a = sum(K(:,:,:,:,1,:,1,1,1,1),6)./sum(K(:,:,:,:,1,:,1,1,1,1)~=0,6);
sens=bart('ecalib -S -m1',a);
figure(111); immontage4D(abs(sens),[])

% sens=ones(size(kspace));
kspace=squeeze(K);

%%
a = sum(sum(kspace,4),5)./sum(sum(kspace~=0,4),5);
sens=bart('ecalib -S -r9 -m1',permute(a,[4 1 2 3]));
sens=sens+1e-7; % no zero vals in sense maps...
immontage4D(abs(sens),[])
%%

sens=sens;
close all;
params=params_init();
params.Lg=2;
params.inspectLg=false
params.L3=5;
params.L4=5;
params.subspacedim2=5;
params.sparsity_transform='TV';
params.Imref=[];
params.x=17;
params.y=50;
params.lambda=2e-2
% sens(sens==0)=1e-2;

params.increase_penalty_parameters=false
params.G.precon=false;
params.G.maxiter=10
P_recon=LRT_recon(kspace,squeeze(sens),params);
%% visualize recon
figure(1000); immontage4D(squeeze(abs(P_recon)),[]);
