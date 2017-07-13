% draft reconstruction code 
clear all; close all; clc;

% 1: make data (most settings in other .m file for now)
uf=0.02; % undersampling factor (excluding center)
noiselevel=0;
ncoils=2;
complexsim=1

% sparsity_transform='wavelet'
sparsity_transform='TVcomplex'

L3=4;               %rank of subspace dimension 3
L4=4;               %rank of subspace dimension 4

run create_undersampled_measurement_S.m    
% >>>>>>>>>>>>>>>>>>>>RECON FROM HERE<<<<<<<<<<<<<<<<<<<<<<<<<<<

% 2: estimate subspaces
nav_parameter_dim1 = squeeze(du(ctrcoords,ctrcoords,:,:,1));
nav_estimate_1= subspace_estimator_multicoil(nav_parameter_dim1,L3);

nav_parameter_dim2 = squeeze(du(ctrcoords,ctrcoords,:,1,:));
nav_estimate_2= subspace_estimator_multicoil(nav_parameter_dim2,L4);

tensorsize=size(du);
imagesize=tensorsize; imagesize(3)=1; 
unfoldedIsize=[size(du,1)*size(du,2),size(du,4)*size(du,5)];                %coil combined
unfoldedKsize=[size(du,1)*size(du,2)*size(du,3),size(du,4)*size(du,5)];     %coils separate

% 3: initialize other operators
if strcmp(sparsity_transform,'wavelet')
    Psi=opWavelet2(res,res,'Daubechies') %wavelet operator (uses SPOT toolbox (+ other dependencies maybe?)
elseif strcmp(sparsity_transform,'TV')
    % to do: can be made ~15 times faster with a finite difference operator
%     Psi=opConvolve(res,res,[-1 1],[0 0],'truncated')* opConvolve(res,res,[-1 1]',[0 0],'truncated') %2D TV operator
elseif strcmp(sparsity_transform,'TVcomplex')
%     Psi=TV2op()
    Psi1=opConvolve(res,res,[-1 1],[0 0],'truncated')
    Psi2=opConvolve(res,res,[-1 1]',[0 0],'truncated')
    Psi=[Psi1;Psi2]
end

sens_normalized=bsxfun(@rdivide,sens,sqrt(sum(abs(sens).^2,3))); 
F=MCFopClass;
set_MCFop_Params(F,(sens_normalized),[res,res],[tensorsize(4),tensorsize(5)]);

%4 zero-filled recon
P0=F'*du;             
P1_0=reshape(P0,unfoldedIsize); %1-unfolding of zero filled recon (what is the shape of this matrix?)
du_1=reshape(du,unfoldedKsize);

figure(4); imshow(abs(P0(:,:,1,1,1)),[]); axis off; title('zero filled recon of one frame')
figure(5); imshow(angle(P0(:,:,1,1,1)),[]); axis off; title('phase of zero filled recon of one frame')

%% 
% lowres_phase_estimate=exp(-1i*angle(P0(:,:,1,1,1)));
% Psi=Psi*opDiag(lowres_phase_estimate(:))

%% ALGO 
%initialize parameters
alpha= 2;         %penalty parameter >0
beta=  2;         %penalty parameter >0
alpha=0.01/(mean(abs(du(du~=0))));beta=alpha; %from paper 

lambda=5e-2;        %sparsity parameter
mu=5e-1 ;           %sparsity parameter

Lg=L3*L4;             %rank of spatial dimension
Lg=3;
niter=10;

%initialize matrices
[Phi,G,C,A,B,Y,Z]= init_G0(P1_0,nav_estimate_1,nav_estimate_2,Lg);                    
A=zeros(size(Psi*G));
Y=zeros(size(A));

MSE=[]; 
for iter=1:niter
    MSE=visualize_convergence(iter,MSE,G,C,Phi,squeeze(I),imagesize,61,33);
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
for ii=1:size(I,3)
    J=[];
    for jj=1:size(I,4);
        J=[J,abs(P_recon(:,:,1,ii,jj))];
    end
    Q=[Q;J];
end
imshow(abs(Q),[0 0.1])
clear Q J 
%%