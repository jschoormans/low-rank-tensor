% qMATCH LOAD all DATA 

addpath(genpath('/opt/amc/bart/')); vars;
%%
cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_09_14\2017_09_14\qM_97389')
file{1}='qm_14092017_1826310_7_2_wip_qmatch_invivo_te70V4.raw'
file{2}='qm_14092017_1828245_8_2_wip_qmatch_invivo_te60V4.raw'
file{3}='qm_14092017_1830051_9_2_wip_qmatch_invivo_te50V4.raw'
file{4}='qm_14092017_1831466_10_2_wip_qmatch_invivo_te40V4.raw'
file{5}='qm_14092017_1833282_11_2_wip_qmatch_invivo_te30V4.raw'
file{6}='qm_14092017_1835098_12_2_wip_qmatch_invivo_te20V4.raw'

 for file_iter=1:6
clear MR
MR=MRecon(file{file_iter})
DTI=0;

aver= repmat([0:199],[8 1]); aver=aver(:);
aver=repmat(aver,[36 1]);

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
% MR.SortData;
%ifft in readout direction + select slice 

read_data_ix = (MR.Parameter.Labels.Index.typ == 1);
ky = MR.Parameter.Labels.Index.ky;
ky_read = ky(read_data_ix);
kz = MR.Parameter.Labels.Index.kz;
kz_read = kz(read_data_ix);
eightline=(mod(aver+4,8)==0);
k0_profile_nr = find((ky_read==0).*(kz_read==0).*eightline);

data = MR.Data(:,k0_profile_nr);
echo = aver(k0_profile_nr);
echo1ch=echo(1:8:end);
echo1ch=(echo1ch+4)/8;
datasorted = zeros(216,25,8);
for i=1:8
    datatemp=data(:,i:8:end);
    for e=1:25
        datasorted(:,echo1ch(e),i)=datatemp(:,e);
    end
end
    

k0_profile{file_iter} = datasorted;

MR.Data=fftshift(ifft(ifftshift(MR.Data,1),[],1),1);
MR.Data=MR.Data(86,:); %132 or 84

aver_tem = MR.Parameter.Labels.Index.aver;
aver_tem(read_data_ix) = aver(:);
MR.Parameter.Labels.Index.aver=aver_tem; clear aver_tem;
K= sortArray(MR);
Kcc=bart('cc -p2',permute(K,[1 2 3 4 7 8 9 10 5 6]));
Kcc=permute(Kcc,[1 2 3 4 9 10 5 6 7 8]);
size(Kcc)
kspace(:,:,:,:,file_iter)=squeeze(Kcc);
size(kspace)
 end
 
kspace=permute(kspace,[6 1 2 3 4 5]);
%%
kspace_binned=reshape(kspace,[213,36,2,8,25,6]); 
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

sens=bart('ecalib -S -r15 -m1',permute(a,[4 1 2 3]));
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
params.L3=5;
params.L4=5;
params.subspacedim1=6;
params.subspacedim2=6; 

Ktemp=cat(4,k0_profile{1},k0_profile{2},k0_profile{3},k0_profile{4},k0_profile{5},k0_profile{6});
size(Ktemp)
Ktemp=permute(Ktemp,[1 5 6 3 2 4]);size(Ktemp)

if 1
    [nav_estimate_1,nav_estimate_2,eigenvals_1,eigenvals_2]= subspace_estimate_3D(Ktemp,params);
    params.nav_estimate_1=nav_estimate_1;
    params.nav_estimate_2=nav_estimate_2;
    params.eigenvals_1=eigenvals_1;
    params.eigenvals_2=eigenvals_2;
end

params.scaleksp=false; 
params.Lg=5;
params.inspectLg=false;
params.sparsity_transform='TV';
params.Imref=[];
params.x=20;
params.y=35;
params.mu=0.02e2;
params.lambda=1.4e-1;
params.automu=0;
params.autolambda=0; 
params.alpha=0.2; 

params.niter=5; 
params.increase_penalty_parameters=false;
params.G.precon=false;
params.G.maxiter=3;

P_recon=LRT_recon(kspace(:,:,:,:,:),squeeze(sens(:,:,:,:)),params);

figure(1000); immontage4D(squeeze(abs(permute(P_recon,[2 1 3 4 5]))),[]);
 
 