clear; close all; clc
if ispc
    cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_08_13')
%     addpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations\tensor\low-rank-tensor'))
else
    cd('/home/qzhang/lood_storage/divi/Ima/parrec/Jasper/LRT/Low_Rank_2017_08_13')
%     addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/tensor/low-rank-tensor'))
%     addpath(genpath('/home/qzhang/lood_storage/divi/Projects/cosart/CS_simulations/tensor/low-rank-tensor'))
    addpath(genpath('/opt/amc/bart/')); vars;
end

MR=MRecon('lr_13082017_1757596_32_2_wip_sc18-vfa-dtiV4.raw')
% MR=MRecon('lr_13082017_1702510_29_2_wip_sc18-vfa-dtiV4.raw')
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
%     MR.Parameter.Parameter2Read.chan=[2,3,4].';
%     MR.Parameter.Parameter2Read.ky=[-50:50].'
    MR.Parameter.Parameter2Read.chan=[34;35;36;37;44;45]

else
    %MR.Parameter.Labels.Index.aver=zeros(size(MR.Parameter.Labels.Index.aver));
    %MR.Parameter.Recon.ImmediateAveraging='Yes'
    MR.Parameter.Recon.ArrayCompression='No';
    MR.Parameter.Recon.ACNrVirtualChannels=6;
    MR.Parameter.Parameter2Read.typ = 1;
    MR.Parameter.Parameter2Read.chan=[34;35;36;37;44;45]
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


[nav_estimate_1,eigenvals_1,nav_estimate_2,eigenvals_2] = full_subspace_estimate_3D(K,params)
%}

%%


K=MR.Data;
K=fftshift(ifft(ifftshift(K,1),[],1),1); 
size(K)
%%
Ktemp=squeeze(MR.Data); 

%%
kspace=K(32,2:end,2:end,:,:,:,:,:,:,:,:,:);
imshow(squeeze(abs(kspace(1,:,:,1,1))),[0 1e-2])
size(kspace)


% remove stupid checkerboard pattern
che=create_checkerboard([1,size(kspace,2),size(kspace,3)]);
kspace=bsxfun(@times,kspace,che);
kspace=squeeze(kspace);

if DTI
par_dim1 = 9; par_dim2 = 6;
du = zeros(size(kspace,1),size(kspace,2),size(kspace,3),par_dim1,par_dim2);
du_temp = zeros(size(Ktemp,1),size(Ktemp,2),size(Ktemp,3),size(Ktemp,4),par_dim1,par_dim2);
for ii = 1:par_dim1
    for jj = 1:par_dim2
        du(:,:,:,ii,jj) = kspace(:,:,:,((jj-1) * par_dim1 +ii));
        du_temp(:,:,:,:,ii,jj) = Ktemp(:,:,:,:,((jj-1) * par_dim1 +ii));
    end
end
kspace=du; 
Ktemp = du_temp;
end

a = sum(sum(kspace(:,:,:,:,:),4),5)./sum(sum(kspace(:,:,:,:,:)~=0,4),5);
sens=bart('ecalib -S -r15 -m1',permute(a,[4 1 2 3]));
sens=sens+1e-7; % no zero vals in sense maps...

figure(112)
subplot(311)
immontage4D(abs(sens),[])
subplot(312)
immontage4D(real(sens),[])
subplot(313)
immontage4D(angle(sens),[-pi pi])

%%

% sens=ones(size(sens)); % WERKT BETER DAN SENSE MAPS 
% sens=abs(sens); % VALT OOK TE PROBEREN 
sens_onecoil=sens(:,:,:,[3]); 
kspace_onecoil=kspace(:,:,[3],:,:);

params=params_init();
params.L3=3;
params.L4=3;
params.subspacedim1=1;
params.subspacedim2=6; 
[nav_estimate_1,nav_estimate_2,eigenvals_1,eigenvals_2]= subspace_estimate_3D(Ktemp(:,:,:,:,:,:),params);

params.nav_estimate_1=nav_estimate_1;
params.nav_estimate_2=nav_estimate_2;
params.eigenvals_1=eigenvals_1;
params.eigenvals_2=eigenvals_2;

params.Lg=3;
params.inspectLg=false;
params.sparsity_transform='TV';
params.Imref=[];
params.x=50;
params.y=35;
params.mu=0.1e2;
params.lambda=5e-3;
% sens(sens==0)=1e-2;

params.niter=15; 
params.increase_penalty_parameters=false;
params.G.precon=true;
params.G.maxiter=50;

P_recon=LRT_recon(kspace_onecoil,squeeze(sens_onecoil),params);
% P_recon=LRT_recon(kspace,squeeze(sens),params);

%% visualize recon
figure(1000); immontage4D(squeeze(abs(P_recon)),[]);

