% recon 28 july
clear all; close all; clc
cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_07_28')  
%%
% MR=MRecon('lr_28072017_1118286_15_2_wip_vfa-t2prep_tenr18dynsV4.raw'); % wrong orientation/resolution/t2prepvals
% MR1=MRecon('lr_28072017_1242319_23_2_wip_vfa-t2prep_tenr18dynsV4.raw')
MR1=MRecon('lr_28072017_1150231_18_2_wip_vfa-dtiV4.raw')
MR1.Perform
%%
% MR=MRecon('lr_28072017_1141160_17_2_wip_vfa-t2prep_tenr18dynsV4.raw')
% MR=MRecon('lr_28072017_1112088_14_2_wip_vfa-t2prep_tenr18dynsV4.raw')
% MR=MRecon('lr_28072017_1242319_23_2_wip_vfa-t2prep_tenr18dynsV4.raw')
MR=MRecon('lr_28072017_1150231_18_2_wip_vfa-dtiV4.raw')

DTI=true;
if ~DTI
MR.Parameter.Labels.Index.aver=(MR.Parameter.Labels.Index.rf);
MR.Parameter.Parameter2Read.aver=[0:60].'
else
    MR.Parameter.Labels.Index.aver=zeros(size(MR.Parameter.Labels.Index.aver));
end
% MR.Parameter.Parameter2Read.aver=[0:60].';
% MR.Parameter.Labels.NumberOfEchoes=60
MR.Parameter.Recon.RemoveMOversampling='No'
MR.Parameter.Recon.RemovePOversampling='No'
MR.Parameter.Recon.ArrayCompression='No'
MR.Parameter.Recon.ACNrVirtualChannels=4
MR.Parameter.Parameter2Read.typ = 1;
MR.Parameter.Recon.ImmediateAveraging='No'

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
K=ifftshift(ifft(K,[],1),1); size(K)
K=K(65,:,:,:,:,:,:,:,:,:,:,:);
% remove stupid checkerboard pattern
che=create_checkerboard([1,size(K,2),size(K,3)]);
K=bsxfun(@times,K,che);
a = sum(K(:,:,:,:,1,:,1,1,1,1),6)./sum(K(:,:,:,:,1,:,1,1,1,1)~=0,6);
sens=bart('ecalib -S -m1',a);
% sens=ones(size(kspace));
%%
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
close all;
params=params_init();
params.Lg=6;
params.inspectLg=true
params.L3=4;
params.L4=4;
params.sparsity_transform='TV';
params.Imref=[];
params.x=30;
params.y=30;

params.increase_penalty_parameters=true
P_recon=LRT_recon(kspace,squeeze(sens),params);