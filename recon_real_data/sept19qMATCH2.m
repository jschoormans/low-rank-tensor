% qMATCH LOAD all DATA 
clear all
addpath(genpath('/opt/amc/bart/')); vars;
%%
% cd('L:\basic\divi\Ima\parrec\Jasper\LRT\2017_09_19\qM_101945')
% file='qm_19092017_2052467_15_2_wip_qmatch_invivo_allteV4.raw'

cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_09_21')
file='lr_21092017_1756316_4_2_wipqmatchinvivote70V4.raw'


clear MR
MR=MRecon(file)
DTI=0;
aver= repmat([0:199],[13 1]); aver=aver(:);
aver=repmat(aver,[36 1]);
aver=repmat(aver,[6 1]); 

% MR.Parameter.Parameter2Read.chan=[2,3].'
MR.Parameter.Parameter2Read.typ = 1;
MR.Parameter.Recon.ImmediateAveraging='No';

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

read_data_ix = (MR.Parameter.Labels.Index.typ == 1);    
% %%%
% ky = MR.Parameter.Labels.Index.ky;
% ky_read = ky(read_data_ix);
% kz = MR.Parameter.Labels.Index.kz;
% kz_read = kz(read_data_ix);
% dyn=MR.Parameter.Labels.Index.dyn;
% dyn_read=dyn(read_data_ix);
% eightline=(mod(aver,8)==3);
% k0_profile_nr = find((abs(ky_read)<3).*(abs(kz_read)<3).*eightline);
% % 
% data = MR.Data(:,k0_profile_nr);
% echo = aver(k0_profile_nr);
% echo1ch=echo(1:8:end);
% echo1ch=(echo1ch+5)/8;
% 
% dyn_read_profile0=dyn_read(k0_profile_nr);
% 
% datasorted = zeros(216,25,6,8);
% for q=1:max(dyn+1);
%     datatemp=data(:,dyn_read_profile0==(q-1));
%     for i=1:8
%         datatemp=datatemp(:,i:8:end);
%         for e=1:25
%             datasorted(:,ky,kz,echo1ch(e),q,i)=datatemp(:,e);
%         end
%     end
% end
% 
% %%%


MR.Data=fftshift(ifft(ifftshift(MR.Data,1),[],1),1);
MR.Data=MR.Data(140,:); %132 or 84

aver_tem = MR.Parameter.Labels.Index.aver;
aver_tem(read_data_ix) = aver(:);

MR.Parameter.Labels.Index.aver=aver_tem; clear aver_tem;

MR.SortData
size(MR.Data)
K=MR.Data;
K=squeeze(K);
% K= sortArray(MR);
% Kcc=bart('cc -p4',permute(K,[1 2 3 4 7 8 9 10 5 6]));
% Kcc=permute(Kcc,[1 2 3 4 9 10 5 6 7 8]);
% size(Kcc)
%%
kspace=squeeze(K);
size(kspace)
kspace=permute(kspace,[6 1 2 3 5 4]);
size(kspace)
kspace_binned=reshape(kspace,[240,37,13,8,25,6]); 
kspace_binned=mean(kspace_binned,4)./(eps+sum(abs(kspace_binned)~=0,4)); 
kspace_binned=squeeze(kspace_binned);
size(kspace_binned)
kspace=permute(kspace_binned,[6 1 2 3 4 5]);

 %%
imshow(squeeze(abs(kspace(1,:,:,1,1))),[0 1e-2])
size(kspace)

% remove stupid checkerboard pattern
che=create_checkerboard([1,size(kspace,2),size(kspace,3)]);
kspace=bsxfun(@times,kspace,che);
kspace=squeeze(kspace);

a = sum(sum(kspace(:,:,:,:,:),4),5)./(sum(sum(kspace(:,:,:,:,:)~=0,4),5)+eps);
a2 = (sum(kspace(:,:,:,:,:),5))./((sum(kspace(:,:,:,:,:)~=0,5))+eps);

% sens=bart('ecalib -S -m1',permute(a,[4 1 2 3]));
sens=bart('caldir 10',permute(a,[4 1 2 3]));

% sens=sens+1e-7; % no zero vals in sense maps...

figure(112)
subplot(311)
immontage4D(abs(sens),[])
subplot(312)
immontage4D(real(sens),[])
subplot(313)
immontage4D(angle(sens),[-pi pi])
%% pre-estimate subspaces


%%
params=params_init();
params.L3=3;
params.L4=4;
params.subspacedim1=6;
params.subspacedim2=3;


if 0
    % Ktemp=cat(4,k0_profile{1},k0_profile{2},k0_profile{3},k0_profile{4},k0_profile{5},k0_profile{6});
    % size(Ktemp)
    % Ktemp=permute(Ktemp,[1 5 6 3 2 4]);size(Ktemp)
    [nav_estimate_1,nav_estimate_2,eigenvals_1,eigenvals_2]= subspace_estimate_3D(Ktemp,params);
    params.nav_estimate_1=nav_estimate_1;
    params.nav_estimate_2=nav_estimate_2;
    params.eigenvals_1=eigenvals_1;
    params.eigenvals_2=eigenvals_2;
end

params.scaleksp=false;
params.Lg=2;
params.inspectLg=false;
params.sparsity_transform='TV';
params.Imref=[];
params.x=20;
params.y=105;
params.mu=0.02e2;
params.lambda=1.4e-2;
params.automu=1;
params.autolambda=1; 
params.alpha=5; 

params.niter=12; 
params.increase_penalty_parameters=false;
params.G.precon=false;
params.G.maxiter=50;
params.Lg=10
coil=2
P_recon=LRT_recon_test(kspace(:,:,:,:,:),squeeze(sens(:,:,:,:)),params);

figure(1000); immontage4D(squeeze(abs(permute(P_recon,[2 1 3 4 5]))),[]);

