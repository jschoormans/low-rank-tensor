% draft reconstruction code 
clear all; close all; clc;

% 1: make data (settings in other .m file for now)
uf=0.2; % undersampling factor (excluding center)
noiselevel=0.2;

run create_undersampled_measurement.m
% 1a: add noise
if noiselevel>0
du=du+(randn(size(du)).*mean(du(:)).*noiselevel).*(du~=0);
end

% 2: estimate subspaces
L3=4;               %rank of subspace dimension 3
L4=4;               %rank of subspace dimension 4

nav_parameter_dim1 = squeeze(du(ctrcoords,ctrcoords,:,5));
nav_estimate_1= subspace_estimator(nav_parameter_dim1,L3);

nav_parameter_dim2 = squeeze(du(ctrcoords,ctrcoords,5,:));
nav_estimate_2= subspace_estimator(nav_parameter_dim2,L4);

% 3: initialize other operators
tensorsize=size(du)
unfoldedsize=[size(du,1)*size(du,2),size(du,3)*size(du,4)]

F=Fop([res,res]);
Psi=opWavelet2(res,res,'Daubechies') %wavelet operator (uses SPOT toolbox (+ other dependencies maybe?) 

%4 zero-filled recon
P0=F'*du;
P1_0=reshape(P0,unfoldedsize); %1-unfolding of zero filled recon (what is the shape of this matrix?)
du_1=reshape(du,unfoldedsize);

%% ALGO 
%initialize parameters
alpha= 1e1;         %penalty parameter >0
beta=  0.1;         %penalty parameter >0
lambda=1;        %sparsity parameter
mu=1e-3 ;            %sparsity parameter
Lg=101;             %rank of spatial dimension
niter=5;

%initialize matrices
[Phi,G,C,A,B,Y,Z]= init_G0(P1_0,nav_estimate_1,nav_estimate_2,Lg);                    

MSE=[]; 
for iter=1:niter
    MSE=visualize_convergence(iter,MSE,G,C,Phi,I,tensorsize,70,70)
    
    Ak=soft_thresh_A(G,Y,alpha,lambda,Psi);             %15
    Bk=soft_thresh_B(C,Z,mu,beta);                      %16
    Gk=conj_grad_G_2(G,C,Ak,Y,alpha,Psi,du_1,Phi,F);     %17
    Ck=C; %temp
%     Ck=conj_grad_C_2(Gk,C,Bk,Z,beta,du_1,Phi,F); % to do...
    Yk=Y+alpha*(Ak-Psi*Gk);
    Zk=Z-beta.*(Bk-Ck);
    
    A=Ak; B=Bk; G=Gk; C=Ck; Y=Yk; Z=Zk; %update iteration
end


%%