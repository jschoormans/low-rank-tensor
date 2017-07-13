clear all; close all; clc
cd('/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/Low_Rank_2017_07_11')
addpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations\tensor\low-rank-tensor'))
%%
% for first experiments we ignore the coil dimension, and do reconstruction
% on a coil-combined fft of an image;, instead of on the k-space directly 
MR_data = MRecon('lo_12072017_1926046_17_2_wipvfat2preplowrankV4.raw');
MR_data.Parameter.Recon.RemoveMOversampling='No'
MR_data.Parameter.Recon.RemovePOversampling='No'
MR_data.Parameter.Parameter2Read.typ = 1;
MR_data.Parameter.Parameter2Read.Update;
% load data
MR_data.ReadData;
MR_data.RandomPhaseCorrection;
MR_data.RemoveOversampling;
MR_data.PDACorrection;
MR_data.DcOffsetCorrection;
MR_data.MeasPhaseCorrection;

MR_data.Parameter.Labels.Index.card=MR_data.Parameter.Labels.Index.rf; %use cardiac dimensions for TFE number
MR_data.SortData;
% TO DO: SORT SUCH THAT ECHO NUMBER IS A DIMENSION
% use MR_data.Parameter.Labels.Index.rf
MR_data.K2IM;

%% 

K=MR_data.Data(128,:,:,:,:,:,:,:,:,:,:,:);
size(K)

% remove stupid checkerboard pattern
che=create_checkerboard([1,132,65]);
K=bsxfun(@times,K,che);

a = sum(K(:,:,:,:,1,:,1,1,1,1),6)./sum(K(:,:,:,:,1,:,1,1,1,1)~=0,6);
sens=bart('ecalib -r 20 -m1',a);
du=squeeze(K);
mask=squeeze(du(:,:,1,:,:))~=0;


figure(1); Q=[];
for ii=1:size(mask,3)
    J=[];
    for jj=1:size(mask,4);
        J=[J,abs(mask(:,:,ii,jj))];
    end
    Q=[Q;J];
end
figure(1); clf
imshow(abs(Q),[0 1])
xlabel('Prameter 1')
ylabel('Prameter 2')
clear Q J 
%% estimate subspaces
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/tensor/low-rank-tensor'))

L3=4;
L4=3;

ctr=5
ctrcoords1=61-ctr:61+ctr %???
ctrcoords2=32-ctr:32+ctr
% >>>>>>>>>>>>>>>>>>>>RECON FROM HERE<<<<<<<<<<<<<<<<<<<<<<<<<<<

% 2: estimate subspaces
% nav_parameter_dim1 = squeeze(du(ctrcoords,ctrcoords,:,:,1));
nav_parameter_dim1 = squeeze(du(ctrcoords1,ctrcoords2,:,:,1));
nav_estimate_1= subspace_estimator_multicoil(nav_parameter_dim1,L3);

nav_parameter_dim2 = squeeze(du(ctrcoords1,ctrcoords2,:,1,:));
nav_parameter_dim2=squeeze(nav_parameter_dim2);
% 

nav_estimate_2= subspace_estimator_multicoil(nav_parameter_dim2,L4);
figure(2); clf
plot(real(nav_estimate_1)); hold on
% plot(imag(nav_estimate_1),'k'); hold off
figure(3); clf
plot(real(nav_estimate_2),'r'); hold on
plot(imag(nav_estimate_2),'k'); hold off

tensorsize=size(du);
imagesize=tensorsize; imagesize(3)=1; 
unfoldedIsize=[size(du,1)*size(du,2),size(du,4)*size(du,5)];                %coil combined
unfoldedKsize=[size(du,1)*size(du,2)*size(du,3),size(du,4)*size(du,5)];     %coils separate

sens=squeeze(sens)
res1=132
res2=65
sparsity_transform='wavelet'
% 3: initialize other operators
if strcmp(sparsity_transform,'wavelet')
    Psi=opWavelet2(res1,res2,'Daubechies') %wavelet operator (uses SPOT toolbox (+ other dependencies maybe?)
elseif strcmp(sparsity_transform,'TV')
    % to do: can be made ~15 times faster with a finite difference operator
    Psi=opConvolve(res,res,[-1 1],[0 0],'truncated')* opConvolve(res,res,[-1 1]',[0 0],'truncated') %2D TV operator
end

sens_normalized=bsxfun(@rdivide,sens,sqrt(sum(abs(sens+eps).^2,3))); 
F=MCFop([res1,res2],(sens_normalized));

%4 zero-filled recon
P0=F'*du;             
P1_0=reshape(P0,unfoldedIsize); %1-unfolding of zero filled recon (what is the shape of this matrix?)
du_1=reshape(du,unfoldedKsize);

figure(4); imshow(abs(P0(:,:,1,2,2))./max(abs(P0(:))),[]); axis off; title('zero filled recon of one frame')
%% ALGO 
%initialize parameters
alpha= 2;         %penalty parameter >0
beta=  2;         %penalty parameter >0
% alpha=0.01/(mean(abs(du(du~=0))));beta=alpha; %from paper 

lambda=5e-1;        %sparsity parameter
mu=5e-1 ;           %sparsity parameter

Lg=L3*L4;             %rank of spatial dimension
Lg=5;
niter=10;

%initialize matrices
[Phi,G,C,A,B,Y,Z]= init_G0(P1_0,nav_estimate_1,nav_estimate_2,Lg);                    
Y=zeros(size(Psi*G)); % to do: fix in init 
MSE=[]; 
for iter=1:niter
    MSE=visualize_convergence(iter,MSE,G,C,Phi,[],imagesize,17,39);
    Ak=soft_thresh_A(G,Y,alpha,lambda,Psi);                     %15
    Bk=soft_thresh_B(C,Z,mu,beta);                              %16
    Gk=precon_conj_grad_G(G,C,Ak,Y,alpha,Psi,du_1,Phi,F);       %17
    Ck=precon_conj_grad_C(Gk,C,Bk,Z,beta,du_1,Phi,F);           %18
    Yk=Y+alpha*(Ak-Psi*Gk);
    Zk=Z+beta.*(Bk-Ck);
    
    G=Gk; C=Ck; Y=Yk; Z=Zk; %update iteration
    
% alpha=alpha*1.5; beta=beta*1.5;    %optional: increasing beta & alpha
end
%%
P_recon=G*C*Phi;

P_recon=reshape(P_recon,imagesize);

figure(998); 
 Q=[]
for ii=1:size(P_recon,4)
    J=[];
    for jj=1:size(P_recon,5);
        J=[J,angle(P_recon(:,:,1,ii,jj))];
    end
    Q=[Q;J];
end
figure(100)
imshow((Q),[])
clear Q J 

%% 

figure(998); 
 Q=[]
for ii=15:24
    J=[];
    for jj=1:size(P_recon,5);
        J=[J,angle(P0(:,:,1,ii,jj))];
    end
    Q=[Q;J];
end
figure(101)
imshow((Q),[])
clear Q J 





