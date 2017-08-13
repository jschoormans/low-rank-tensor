% recon 28 july
clear all; close all; clc
cd('/home/qzhang/lood_storage/divi/Ima/parrec/Kerry/Data/2017_07_31_LRT_DTI_T2prep');
%
% MR=MRecon('lr_28072017_1118286_15_2_wip_vfa-t2prep_tenr18dynsV4.raw'); % wrong orientation/resolution/t2prepvals
% MR1=MRecon('lr_28072017_1242319_23_2_wip_vfa-t2prep_tenr18dynsV4.raw')
% MR1=MRecon('lr_28072017_1150231_18_2_wip_vfa-dtiV4.raw')
% MR1.Perform

% MR=MRecon('lr_28072017_1141160_17_2_wip_vfa-t2prep_tenr18dynsV4.raw')
% MR=MRecon('lr_28072017_1112088_14_2_wip_vfa-t2prep_tenr18dynsV4.raw')
% MR=MRecon('lr_28072017_1242319_23_2_wip_vfa-t2prep_tenr18dynsV4.raw')
MR=MRecon('lr_31072017_1737345_2_2_wip_sc18-vfa-dtiV4.raw')

DTI=true;

if ~DTI
    MR.Parameter.Labels.Index.aver=(MR.Parameter.Labels.Index.rf);
    MR.Parameter.Parameter2Read.aver=[0:max(MR.Parameter.Labels.Index.aver)].'
    % MR.Parameter.Parameter2Read.aver=[0:60].';
    % MR.Parameter.Labels.NumberOfEchoes=60
    % MR.Parameter.Recon.RemoveMOversampling='No'
    % MR.Parameter.Recon.RemovePOversampling='No'
    MR.Parameter.Recon.ArrayCompression='No';
    MR.Parameter.Recon.ACNrVirtualChannels=2;
    MR.Parameter.Parameter2Read.typ = 1;
    MR.Parameter.Parameter2Read.ky=[-28:29].'
    MR.Parameter.Parameter2Read.kz=[-28:29].'

MR.Parameter.Recon.ImmediateAveraging='No';
MR.Parameter.Parameter2Read.chan=[2,3].';
else
%     MR.Parameter.Labels.Index.aver=zeros(size(MR.Parameter.Labels.Index.aver));
%     MR.Parameter.Recon.ImmediateAveraging='Yes'
MR.Parameter.Labels.Index.aver = 0 .* MR.Parameter.Labels.Index.rf;
MR.Parameter.Parameter2Read.typ = 1;
% run retro_profile_check; %optional

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

K=MR.Data;
K=ifftshift(ifft(K,[],1),1); size(K)
K=K(32,:,:,:,:,:,:,:,:,:,:,:);
% remove stupid checkerboard pattern
che=create_checkerboard([1,size(K,2),size(K,3)]);
K=bsxfun(@times,K,che);
a = sum(K(:,:,:,:,1,:,1,1,1,1),6)./sum(K(:,:,:,:,1,:,1,1,1,1)~=0,6);
% sens=ones(size(kspace));

kspace=squeeze(K);
if ~DTI
if ndims(kspace)==4 % if one chan
   kspace=permute(kspace,[1 2 5 3 4]); 
end
end

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
%% crop kspace
kspace=du; 
kspace=bart('cc -p4',permute(kspace,[6 1 2 3 7 4 5]));
kspace=permute(kspace,[2 3 4 6 7 1 5]);
kspace=kspace(3:end-2,2:end-1,:,:,:);
% a=sum(sum(kspace,4),5)./sum(sum(kspace~=0,4),5);
a=kspace(:,:,:,1,1);
sens=bart('ecalib -S -m1',permute(a,[4 1 2 3]));
figure(1002); immontage4D(abs(squeeze(sens)),[])

%%
close all;
params=params_init();
params.Lg=1;
params.inspectLg=false;
params.L3=5;
params.L4=5;
params.sparsity_transform='TV';
params.Imref=[];
params.x=25;
params.y=25;
params.G.maxiter=10;
params.lambda=4e-3;
params.mu=1e3;

params.increase_penalty_parameters=false
P_recon=LRT_recon(kspace,squeeze(sens),params);

%%
figure(1000); immontage4D(abs(squeeze(P_recon)),[])
figure(1001); immontage4D(angle(squeeze(P_recon)),[])
figure(1002); surf(squeeze(mean(mean(abs(squeeze(P_recon(20:30,20:30,1,:,:))),1),2)))