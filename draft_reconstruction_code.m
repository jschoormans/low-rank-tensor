% draft reconstruction code 
clear all; close all; clc;

% 1: make data (most settings in other .m file for now)

uf=0.01; % undersampling factor (excluding center)
noiselevel=0;
sparsity_transform='wavelet'
% sparsity_transform='TV'

run create_undersampled_measurement.m    


% 2: estimate subspaces
L3=4;               %rank of subspace dimension 3
L4=6;               %rank of subspace dimension 4

nav_parameter_dim1 = squeeze(du(ctrcoords,ctrcoords,:,1));
nav_estimate_1= subspace_estimator(nav_parameter_dim1,L3);

nav_parameter_dim2 = squeeze(du(ctrcoords,ctrcoords,1,:));
nav_estimate_2= subspace_estimator(nav_parameter_dim2,L4);

% 3: initialize other operators
tensorsize=size(du);
unfoldedsize=[size(du,1)*size(du,2),size(du,3)*size(du,4)];

F=Fop([res,res]);

if strcmp(sparsity_transform,'wavelet')==1
Psi=opWavelet2(res,res,'Daubechies') %wavelet operator (uses SPOT toolbox (+ other dependencies maybe?) 
else
Psi=opConvolve(res,res,[-1 1],[0 0],'truncated')* opConvolve(res,res,[-1 1]',[0 0],'truncated') %2D TV operator
end

%4 zero-filled recon
P0=F'*du;
P1_0=reshape(P0,unfoldedsize); %1-unfolding of zero filled recon (what is the shape of this matrix?)
du_1=reshape(du,unfoldedsize);

figure(4); imshow(abs(P0(:,:,2,2)),[])
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
    Gk=conj_grad_G_3(G,C,Ak,Y,alpha,Psi,du_1,Phi,F);    %17 to do...
    Ck=conj_grad_C_3(Gk,C,Bk,Z,beta,du_1,Phi,F);        %18 to do...
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