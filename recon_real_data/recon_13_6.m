clear all; close all; clc
raw_data_fn = 'L:\basic\divi\Ima\parrec\Jasper\Low_Rank_2017_06_13\lo_13062017_1942169_11_2_wipdptset2prepallinonedecendingt2acendingbV4.raw';
MR_data_raw = MRecon(raw_data_fn);

%%
% for first experiments we ignore the coil dimension, and do reconstruction
% on a coil-combined fft of an image;, instead of on the k-space directly 
MR_data = MR_data_raw.Copy;
MR_data.Parameter.Recon.RemoveMOversampling='No'
MR_data.Parameter.Recon.RemovePOversampling='No'
MR_data.Perform
MR_data.ShowData
%%
ima_data_cc = squeeze(double(MR_data.Data));
ima_data_cc_rs = reshape(ima_data_cc, 640, 320, 10, 9);
kspa_fs = bart('fft 3', ima_data_cc_rs(160:479,:,:,:));
res = size(kspa_fs, 1);

%% undersample

%% Simulate undersampling
I=ima_data_cc_rs(160:479,:,:,:);
F=Fop([res,res]);
kspa_fs = F*I;

%---gen mask
uf = 0.2;
mask=rand(size(kspa_fs))>(1-uf); %undersampling

% add center for subspace estimae
ctr=10;
ctrcoords=floor(res/2)-ctr: floor(res/2)+ctr; 
ll=length(ctrcoords);
mask(ctrcoords,ctrcoords,1,:)=ones(ll,ll,1,size(kspa_fs,4));
mask(ctrcoords,ctrcoords,:,1)=ones(ll,ll,size(kspa_fs,3),1);

% for all measurements
ctrcoordsall=floor(res/2)-3: floor(res/2)+3; 
mask(ctrcoordsall,ctrcoordsall,:,:)=ones(7,7,size(kspa_fs,3),size(kspa_fs,4));

figure(3); Q=[]
for ii=1:size(kspa_fs,3)
    J=[];
    for jj=1:size(kspa_fs,4);
        J=[J,abs(mask(:,:,ii,jj))];
    end
    Q=[Q;J];
end
imshow(abs(Q),[0 1]); clear Q J 

kspa_us = kspa_fs.*mask;
%% RECON

sparsity_transform='wavelet'
% sparsity_transform='TV'

% 2: estimate subspaces
L3=4;               %rank of subspace dimension 3
L4=6;               %rank of subspace dimension 4

nav_parameter_dim1 = squeeze(kspa_us(ctrcoords,ctrcoords,:,1));
nav_estimate_1= subspace_estimator(nav_parameter_dim1,L3);

nav_parameter_dim2 = squeeze(kspa_us(ctrcoords,ctrcoords,1,:));
nav_estimate_2= subspace_estimator(nav_parameter_dim2,L4);

% 3: initialize other operators
tensorsize=size(kspa_us);
unfoldedsize=[size(kspa_us,1)*size(kspa_us,2),size(kspa_us,3)*size(kspa_us,4)];


if sparsity_transform=='wavelet'
Psi=opWavelet2(res,res,'Daubechies') %wavelet operator (uses SPOT toolbox (+ other dependencies maybe?) 
else
Psi=opConvolve(res,res,[-1 1],[0 0],'truncated')* opConvolve(res,res,[-1 1]',[0 0],'truncated') %2D TV operator
end

%4 zero-filled recon
P0=F'*kspa_us;
P1_0=reshape(P0,unfoldedsize); %1-unfolding of zero filled recon (what is the shape of this matrix?)
kspa_us_1=reshape(kspa_us,unfoldedsize);

figure(4); imshow(abs(P0(:,:,1,1)),[])
%% ALGO 
%initialize parameters
alpha= 0.2;         %penalty parameter >0
beta=  0.2;         %penalty parameter >0
lambda=5e-2;        %sparsity parameter
mu=1e-2 ;           %sparsity parameter
Lg=24;             %rank of spatial dimension
niter=50;

%initialize matrices
[Phi,G,C,A,B,Y,Z]= init_G0(P1_0,nav_estimate_1,nav_estimate_2,Lg);                    

MSE=[]; 
for iter=1:niter
    MSE=visualize_convergence(iter,MSE,G,C,Phi,I,tensorsize,210,170)
    
    Ak=soft_thresh_A(G,Y,alpha,lambda,Psi);             %15
    Bk=soft_thresh_B(C,Z,mu,beta);                      %16
    Gk=conj_grad_G_3(G,C,Ak,Y,alpha,Psi,kspa_us_1,Phi,F);    %17 to do...
    Ck=conj_grad_C_3(Gk,C,Bk,Z,beta,kspa_us_1,Phi,F);        %18 to do...
    Yk=Y+alpha*(Ak-Psi*Gk);
    Zk=Z+beta.*(Bk-Ck);
    
    G=Gk; C=Ck; Y=Yk; Z=Zk; %update iteration
end
P_recon=G*C*Phi;
P_recon=reshape(P_recon,tensorsize);
%%
figure(998); 
 Q=[]
for ii=1:size(I,3)
    J=[];
    for jj=1:size(I,4);
        J=[J,abs(P_recon(:,:,ii,jj))];
    end
    Q=[Q;J];
end
imshow(abs(Q),[0 0.1])
clear Q J 
%%










