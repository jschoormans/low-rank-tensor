% draft reconstruction code 
clear all; close all; clc;

% 1: make data (settings in other .m file for now)
run create_undersampled_measurement.m

% 2: estimate subspaces
L3=4; %rank of subspace dimension 3
L4=4; %rank of subspace dimension 3

nav_parameter_dim1 = squeeze(du(ctrcoords,ctrcoords,:,5));
nav_estimate_1= subspace_estimator(nav_parameter_dim1,L3);

nav_parameter_dim2 = squeeze(du(ctrcoords,ctrcoords,5,:));
nav_estimate_2= subspace_estimator(nav_parameter_dim2,L4);

% 3: initialize other operators
F=Fop(size(du),[res*res,size(du,3)*size(du,4)]);
Psi=opWavelet2(res,res,'Daubechies') %wavelet operator (uses SPOT toolbox (+ other dependencies maybe?) 

%4 zero-filled recon
P0=F'*du;
P1_0=reshape(P0,[res*res,size(du,3)*size(du,4)]); %1-unfolding of zero filled recon (what is the shape of this matrix?)
%% ALGO 
%initialize parameters
alpha= 0.1;         %penalty parameter >0
beta=  0.1;         %penalty parameter >0
lambda=1e-4;        %sparsity parameter
mu=1e-3             %sparsity parameter
Lg=100;             %rank of spatial dimension
niter=5;

% add noise to test l2 norm gradient (TEMP)
P1_0=P1_0+rand(size(P1_0)).*mean(P1_0(:)).*20;

%initialize matrices
Phi=kron(nav_estimate_2,nav_estimate_1);    %from subspaces
G= init_G0(P1_0,Phi,Lg);                    % first Lg vectors from left-dominant  svd of P10 Psi^H
C=G'*P1_0*(Phi);                        % G0^H P10 Psi^H 

A = zeros(size(G));
B = zeros(size(C));
Y = zeros(size(G));
Z = zeros(size(C));

MSE=[]; 

for iter=1:niter
    MSE=visualize_convergence(iter,MSE,G,C,Phi',I,size(du),70,70)
    
    Ak=soft_thresh_A(G,Y,alpha,lambda,Psi); %15
    Bk=soft_thresh_B(C,Z,mu,beta); %16
    Gk=conj_grad_G_2(G,C,A,Y,alpha,Psi,du,Phi,F);
    Ck=C; % to do....
    Yk=Y+alpha*(Ak-Psi*Gk);
    Zk=Z-beta.*(Bk-Ck);
    
    A=Ak; B=Bk; G=Gk; C=Ck; Y=Yk; Z=Zk; %update iteration
end



%%