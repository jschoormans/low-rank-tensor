clear; close all; clc
if ispc
<<<<<<< HEAD
    cd('L:\basic\divi\Ima\parrec\Kerry\LRT_Data\2017_09_07')
else
    cd(['/home/',getenv('USER'),'/lood_storage/divi/Ima/parrec/Kerry/LRT_Data/2017_09_07']) %data is sorted in yesterday's folder
end

MR=MRecon('lr_07092017_1815258_2_2_wipsc18dtit2prepcsV4.raw')
=======
    cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_08_15\')
else
    cd(['/home/',getenv('USER'),'/lood_storage/divi/Ima/parrec/Kerry/LRT_Data/2017_09_07/knee']) %data is sorted in yesterday's folder
end

MR=MRecon('lr_07092017_2128474_7_2_wiplrtkneedtit22mmV4.raw')
>>>>>>> c3ad55bff67914f8d4782fff73d87136110af38d
%%
DTI=1;

if(~DTI)
    MR.Parameter.Labels.Index.aver=(MR.Parameter.Labels.Index.rf);
    MR.Parameter.Parameter2Read.aver=[0:max(MR.Parameter.Labels.Index.aver)].';
end
MR.Parameter.Recon.ArrayCompression='No';
MR.Parameter.Recon.ACNrVirtualChannels=4;
MR.Parameter.Parameter2Read.typ = 1;



% load data
disp('readdata')
MR.ReadData;
<<<<<<< HEAD
check_k0_intensity(MR,9,6); %optional
=======
>>>>>>> c3ad55bff67914f8d4782fff73d87136110af38d
MR.RandomPhaseCorrection;
disp('corrections...')
MR.RemoveOversampling;
MR.PDACorrection; %???
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;

disp('sortdata')
MR.SortData;


K=MR.Data;
Ktemp=squeeze(MR.Data);
K=squeeze(K);

K=fftshift(ifft(ifftshift(K,1),[],1),1);

size(K)

%%
<<<<<<< HEAD
kspace=K(30,1:end,1:end,:,:,:,:,:,:,:,:,:);
=======
kspace=K(50,1:end,1:end,:,:,:,:,:,:,:,:,:);
>>>>>>> c3ad55bff67914f8d4782fff73d87136110af38d
%coil compression; 5th dimension must be zeros.
kspace = permute(bart('cc -p 6',permute(kspace, [1 2 3 4 6 5])),[1 2 3 4 6 5]);
imshow(squeeze(abs(kspace(1,:,:,1,1))),[0 1e-2])
size(kspace)

% remove stupid checkerboard pattern
che=create_checkerboard([1,size(kspace,2),size(kspace,3)]);
kspace=bsxfun(@times,kspace,che);
kspace=squeeze(kspace);

if DTI
<<<<<<< HEAD
    par_dim1 = 9; par_dim2 = 6;
=======
    par_dim1 = 7; par_dim2 = 5;
>>>>>>> c3ad55bff67914f8d4782fff73d87136110af38d
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
figure(111); immontage4D(squeeze(abs(kspace(:,:,1,:,:))>0),[]);
%%
% a = sum(sum(kspace(:,:,:,:,:),4),5)./sum(sum(kspace(:,:,:,:,:)~=0,4),5);
<<<<<<< HEAD
a=kspace(:,:,:,1,6);  %highest SNR volume
=======
a=kspace(:,:,:,1,5);  %highest SNR volume
>>>>>>> c3ad55bff67914f8d4782fff73d87136110af38d
sens=bart('ecalib -S -m1 -c0.1',permute(a,[4 1 2 3]));
% sens=sens+1e-7; % no zero vals in sense maps...

figure(112)
subplot(311)
immontage4D(abs(sens),[])
subplot(312)
immontage4D(real(sens),[])
subplot(313)
immontage4D(angle(sens),[-pi pi])

%% 3,4, ok, 2 too noisy
coilnr=[3]
sens_onecoil=sens(:,:,:,coilnr);
kspace_onecoil=kspace(:,:,coilnr,:,:);


params=params_init();
params.L3=3;
params.L4=3;
<<<<<<< HEAD
params.subspacedim1=6;
=======
params.subspacedim1=5;
>>>>>>> c3ad55bff67914f8d4782fff73d87136110af38d
params.subspacedim2=1;
% [nav_estimate_1,nav_estimate_2,eigenvals_1,eigenvals_2]= subspace_estimate_3D(Ktemp(:,:,:,2,:,:),params);
%
% params.nav_estimate_1=nav_estimate_1;
% params.nav_estimate_2=nav_estimate_2;
% params.eigenvals_1=eigenvals_1;
% params.eigenvals_2=eigenvals_2;

<<<<<<< HEAD
params.Lg=3;
=======
params.Lg=5;
>>>>>>> c3ad55bff67914f8d4782fff73d87136110af38d
params.inspectLg=false;
params.sparsity_transform='TV';
params.Imref=[];
params.x=40;
<<<<<<< HEAD
params.y=50;

params.mu=1e1;
params.lambda=5e-3;
params.automu=0;
params.autolambda=0;
=======
params.y=40;

params.mu=1e2;
params.lambda=5e-2;
params.automu=1;
params.autolambda=1;
>>>>>>> c3ad55bff67914f8d4782fff73d87136110af38d

params.scaleksp=0;
params.niter=5;
params.increase_penalty_parameters=false;
params.G.precon=false;
params.G.maxiter=40;
params.normalize_sense=1;

<<<<<<< HEAD
P_recon=LRT_recon(kspace,squeeze((sens)),params);
% P_recon=LRT_recon(kspace_onecoil,squeeze((sens_onecoil)),params);
=======
% P_recon=LRT_recon(kspace,squeeze((sens)),params);
P_recon=LRT_recon(kspace_onecoil,squeeze((sens_onecoil)),params);
>>>>>>> c3ad55bff67914f8d4782fff73d87136110af38d



%% visualize recon
figure(1000); immontage4D(squeeze(abs(P_recon)),[]);
figure(1001);imshow(abs(P_recon(:,:,1,1,6)).',[])
figure(1002)
imshow(cat(1,cat(2,...
    abs(P_recon(:,:,1,1,6).'),...
    abs(P_recon(:,:,1,3,6).'),...
    abs(P_recon(:,:,1,9,6).')),...
    cat(2,abs(P_recon(:,:,1,1,6).'),...
    abs(P_recon(:,:,1,3,6).'),...
    abs(P_recon(:,:,1,9,6).'))),[]);

figure(1003);imshow(abs(sqrt(sum(P_recon(:,:,1,3,:).^2,5))).',[])

%% DTI and T2prep data processing
if(1)
    %%
    DTI_T2_Data =  squeeze(P_recon);
    [y, z, DTI_dim, T2prep_dim] = size(DTI_T2_Data);
    DTI_T2_Data = flipdim(DTI_T2_Data,4);
    
<<<<<<< HEAD
    intensity_cut_off = 10;
=======
    intensity_cut_off = 1;
>>>>>>> c3ad55bff67914f8d4782fff73d87136110af38d
    DTI_T2_Data = DTI_T2_Data.*(abs(DTI_T2_Data)>intensity_cut_off);
    figure(60); immontage4D(abs(DTI_T2_Data),[]);
    
    %%
    %T2 fitting
    TEs = [10 20 40 60 80 100];
<<<<<<< HEAD
    for d =1:DTI_dim
        d
        data = DTI_T2_Data(:,:,d,:);
        [T2_mono_all_allslice(:,:,d), T2_mono_all_rgb_allslice(:,:,:,d), rsquare_mono_all_allslice(:,:,d)] = T2_fitting(data, TEs);
=======
    TEs = [10 30 50 70 90];

    for d =1:DTI_dim
        d
        data = DTI_T2_Data(:,:,d,:);
        [T2_mono_all_allslice(:,:,d), T2_mono_all_rgb_allslice(:,:,:,d), rsquare_mono_all_allslice(:,:,d)] = T2_fitting(data, TEs,1);
>>>>>>> c3ad55bff67914f8d4782fff73d87136110af38d
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
<<<<<<< HEAD
    
    for t =1:T2prep_dim
        t
        DTI_data = abs(permute(DTI_T2_Data(:,:,:,t),[1 2 4 3]));
        [MD(:,:,t), FA(:,:,t), eigvec(:,:,t)] = DTI_fitting(DTI_data, g, b);
    end
    
    figure(63);
    montage(permute(MD,[1 2 4 3]),'displayrange',[0 0.003]); colormap hot; colorbar; title('MD')
    figure(64);
    montage(permute(FA,[1 2 4 3]),'displayrange',[0 1]); colorbar; title('FA')
    iptwindowalign(64,'right',63,'left');
=======
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
    montage(permute(FA,[1 2 4 3]),'displayrange',[0 0.5]);  colorbar; title('FA'); colormap hot;
    figure(65);montage(permute(squeeze(eigvec(:,:,:,:,1)),[1 2 4 3]));  colorbar; title('eigenvector #1');
    figure(66);montage(permute(squeeze(eigvec(:,:,:,:,2)),[1 2 4 3]));  colorbar; title('eigenvector #2');
    figure(67);montage(permute(squeeze(eigvec(:,:,:,:,3)),[1 2 4 3]));  colorbar; title('eigenvector #3');

>>>>>>> c3ad55bff67914f8d4782fff73d87136110af38d
end

