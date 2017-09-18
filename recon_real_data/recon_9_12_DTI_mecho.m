clear; close all; clc
if ispc
    cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_08_13')
    cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_09_05\2017_09_05\lr_23248')
%   addpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations\tensor\low-rank-tensor'))
else
    cd(['/home/',getenv('USER'),'/lood_storage/divi/Ima/parrec/Kerry/LRT_Data/2017_09_12'])
end
%%
clear MR
MR=MRecon('lr_12092017_2004355_2_2_wipdtimechocslowresV4.raw')

DTI=0;
MR.Parameter.Labels.Index.aver=(MR.Parameter.Labels.Index.rf);
MR.Parameter.Parameter2Read.aver=[0:max(MR.Parameter.Labels.Index.aver)].'
MR.Parameter.Recon.ArrayCompression='No';
MR.Parameter.Recon.ACNrVirtualChannels=6;
MR.Parameter.Parameter2Read.typ = 1;
MR.Parameter.Recon.ImmediateAveraging='No';
% MR.Parameter.Parameter2Read.chan=[2;3]
% load data
disp('readdata')
MR.ReadData;
MR.RandomPhaseCorrection;
disp('corrections...')
MR.RemoveOversampling;
MR.PDACorrection; %???
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;
%% MRecon sort
disp('sortdata')
MR.SortData;
K=MR.Data;
Ktemp=squeeze(MR.Data);
K=squeeze(K);
K=fftshift(ifft(ifftshift(K,1),[],1),1);
size(K)

recon_x_loc = 46;
K_1x = K(recon_x_loc,:,:,:,:,:);
Kcc=bart('cc -p5',permute(K_1x,[1 2 3 4 7 8 9 10 5 6]));
Kcc=permute(Kcc,[1 2 3 4 9 10 5 6 7 8]);
size(Kcc)
%% self sort: small memeory allocation
disp('sortdata')
%ifft in readout direction + select slice 
MR.Data=fftshift(ifft(ifftshift(MR.Data,1),[],1),1);

recon_x_loc = 32;
MR.Data=MR.Data(recon_x_loc,:);
K= sortArray(MR);
Kcc=bart('cc -p5',permute(K,[1 2 3 4 7 8 9 10 5 6]));
Kcc=permute(Kcc,[1 2 3 4 9 10 5 6 7 8]);
size(Kcc)

%% all channel recon
kspace=Kcc(1,1:end,1:end,:,:,:,:,:,:,:,:,:);
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

%% one channel recon 
temp=squeeze(Kcc(1,1:end,1:end,:,:,:,:,:,:,:,:,:));

kspace = temp(:,:,1,:,:);

% remove stupid checkerboard pattern
che=create_checkerboard([size(kspace,1),size(kspace,2)]);
kspace=bsxfun(@times,kspace,che);
sens = ones(size(kspace,1), size(kspace, 2));
%%
params=params_init();
params.L3=6;
params.L4=3;
params.subspacedim1=1;
params.subspacedim2=1; 
params.scaleksp=false; 

params.Lg=3;
params.inspectLg=false;
params.sparsity_transform='TV';
params.Imref=[];
params.x=26;
params.y=44;
params.mu=0.1e2;
params.lambda=5e-3;
% sens(sens==0)=1e-2;

params.niter=5; 
params.increase_penalty_parameters=false;
params.G.precon=true;
params.G.maxiter=10;

P_recon=LRT_recon(kspace,squeeze(sens),params);
%% visualize recon
figure(1000); immontage4D(squeeze(abs(P_recon)),[]);
%% DTI and T2prep data processing
if(1)
    %%
    DTI_T2_Data =  squeeze(P_recon);
    [y, z, DTI_dim, T2prep_dim] = size(DTI_T2_Data);
    
    intensity_cut_off = 1;
    DTI_T2_Data = DTI_T2_Data.*(abs(DTI_T2_Data)>intensity_cut_off);
    figure(60); immontage4D(abs(DTI_T2_Data),[]);
    
    %%
    %T2 fitting
    TEs = [0:19]*6.1;
%     TEs = [10 30 50 70 90];

    for d =1:DTI_dim
        d
        data = DTI_T2_Data(:,:,d,:);
        [T2_mono_all_allslice(:,:,d), T2_mono_all_rgb_allslice(:,:,:,d), rsquare_mono_all_allslice(:,:,d)] = T2_fitting(data, TEs,20);
    end
    
    
    figure(61);
    montage(permute(T2_mono_all_allslice,[1 2 4 3]),'displayrange',[]); colormap hot; colorbar
    figure(62);
    montage(permute(rsquare_mono_all_allslice,[1 2 4 3]),'displayrange',[0.9 1]); colormap hot; colorbar
    iptwindowalign(61,'right',62,'left');
    %%
    %DTI
    b = 600;
    g = [0 0 0;...
        0.00000,  0.00000,  1.00000;...
        0.78225,  0.13995,  0.60704;...
        0.60142, -0.66611,  0.44112;...
        -0.29250,  0.74723,  0.59674;...
        -0.11759, -0.90381,  0.41147;...
        -0.92889,  0.28450,  0.23710;...
        0.51744,  0.80103,  0.30102;...
        -0.70473, -0.35950,  0.61166];
    clear MD FA eigvec
    for t =1:T2prep_dim
        t
        DTI_data = abs(permute(DTI_T2_Data(:,:,:,t),[1 2 4 3]));
        [MD(:,:,t), FA(:,:,t), eigvec(:,:,:,:,:,t)] = DTI_fitting(DTI_data, g, b);
    end
    eigvec = permute(eigvec, [1 2 6 4 5 3]);
    
    mask = DTI_data(:,:,:,1)>1;
    MD = bsxfun(@times,MD, mask );
    FA = bsxfun(@times,FA, mask );
    eigvec = bsxfun(@times,eigvec, mask );
    
    figure(63);
    montage(permute(MD,[1 2 4 3]),'displayrange',[]); colormap jet; colorbar; title('MD');
    figure(64);
    montage(permute(FA,[1 2 4 3]),'displayrange',[0 0.7]);  colorbar; title('FA'); colormap hot;
    figure(65);montage(permute(squeeze(eigvec(:,:,:,:,1)),[1 2 4 3]));  colorbar; title('eigenvector #1');
    figure(66);montage(permute(squeeze(eigvec(:,:,:,:,2)),[1 2 4 3]));  colorbar; title('eigenvector #2');
    figure(67);montage(permute(squeeze(eigvec(:,:,:,:,3)),[1 2 4 3]));  colorbar; title('eigenvector #3');

end
