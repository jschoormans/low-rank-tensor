%2D phantom
 clc; close all;

ima_data_cc_rs = reshape(ima_data_cc, 640, 320, 10, 9);
T2Prep = [100 80 70 60 50 40 30 20 10];
% b_value = 900:-100:0;
b_value = 0:50:450;

I = abs(ima_data_cc_rs(160:479,:,:,:));
desired_res = [128 128];

K_full = bart('fft 3', I);
K_downsample = K_full(floor(size(K_full,1)/2-desired_res(1)/2):floor(size(K_full,1)/2-desired_res(1)/2)+desired_res(1)-1,...
    floor(size(K_full,2)/2-desired_res(2)/2):floor(size(K_full,2)/2-desired_res(2)/2)+desired_res(1)-1,:,:);
I = bart('fft -i 3', K_downsample);
clear K_full desired_res K_downsample

I = abs(I); %use real data for now
%% TUCKER MODEL 
% aa=5:3:65
%  iter=1:length(aa)

rank1=32;
rank2=rank1
rank3=2;
rank4=2;

[F1,F2,F3,F4]=tucker(I,[rank1,rank2,rank3,rank4]);

disp('percentage of original information needed:')
frac=(rank1*rank2*rank3*rank4)./prod(size(I))
disp('frct of orignal signal explained:' ); 
pct=F3./100


% figure(2)
% plot(frac,pct);
% xlabel('fraction of original signal size')
% ylabel('explained variance by LRT model')

%Check LRT recon
C = F2;
Psi=kron(F1{4},F1{3}); %kronecker prokspa_usct of G^3 and G^4 (echo times and b-values)
G=kron(F1{2},F1{1});  % spatial low-rank matrix 
C1=reshape(C,[rank1*rank2,rank3*rank4]); %flattened C 

image_recon=(G*C1*Psi.');
image_tensor_res=reshape(image_recon,size(I));

image_tensor_res_normalized = image_tensor_res./max(image_tensor_res(:));
I_normalized = I./max(I(:));
error = abs(image_tensor_res_normalized-I_normalized)./(I_normalized+1e-20);
size(error)
figure(7); montage(error(:,:,1,:),'displayrange',[0 0.1]); colormap hot
figure(8); montage(I(:,:,1,:),'displayrange',[]); colormap hot
figure(9); montage(image_tensor_res(:,:,1,:),'displayrange',[]); colormap hot


figure(10); plot(T2Prep,squeeze(image_tensor_res(180,160,10,:)),'b')
hold on; plot(T2Prep,squeeze(I(180,160,10,:)),'r')
figure(11); imagesc(F1{3}); axis off; plot(abs(F1{3}))

%% Simulate undersampling
kspa_fs = bart('fft 3', I);
res = size(kspa_fs, 1);


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
for ii=1:size(I,3)
    J=[];
    for jj=1:size(I,4);
        J=[J,abs(mask(:,:,ii,jj))];
    end
    Q=[Q;J];
end
imshow(abs(Q),[0 1]); clear Q J 

kspa_us = kspa_fs.*mask;

%% RECON


% 1: make data (most settings in other .m file for now)

uf=0.03; % undersampling factor (excluding center)
noiselevel=0;
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

F=Fop([res,res]);

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
Lg=400;             %rank of spatial dimension
niter=50;

%initialize matrices
[Phi,G,C,A,B,Y,Z]= init_G0(P1_0,nav_estimate_1,nav_estimate_2,Lg);                    

MSE=[]; 
for iter=1:niter
    MSE=visualize_convergence(iter,MSE,G,C,Phi,I,tensorsize,61,33)
    
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

