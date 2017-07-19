function LRT_recon(kspace,sens,params,varargin)
% 4D Low-Rank Tensor reconstruction 
% inputs:
% kspace: undersampled kspace with dimensions (x,y,coilsparam1,param2)
% sens  : sens maps - size (x,y,coils) 
% params:;
% optional: I= 'gold truth' tensor image 

% based on IEEE-TMI paper by Jingfei He  ,2016 

% 2017 J Schoormans and Q Zhang - AMC Amsterdam


% Define parameters
sparsity_transform='TVcomplex'
L3=2;
L4=2;

%initialize parameters
alpha= 2;         %penalty parameter >0
beta=  2;         %penalty parameter >0
alpha=0.01/(mean(abs(kspace(kspace~=0))));
beta=alpha; %from paper 

lambda=1e-8;        %sparsity parameter
mu=1e-1;           %sparsity parameter

% Lg=L3*L4;             %rank of spatial dimension
Lg=3;
niter=10;

res1=size(kspace,1);
res2=size(kspace,2);
tensorsize=size(kspace);
imagesize=tensorsize; imagesize(3)=1; 
unfoldedIsize=[size(kspace,1)*size(kspace,2),size(kspace,4)*size(kspace,5)];                %coil combined
unfoldedKsize=[size(kspace,1)*size(kspace,2)*size(kspace,3),size(kspace,4)*size(kspace,5)];     %coils separate

%%
% FIND MASK
mask=squeeze(kspace(:,:,1,:,:))~=0;
figure(1); clf; immontage4D(mask,[0 1]);
xlabel('Parameter 1'); ylabel('Parameter 2');


% TO DO: FIND COORDINATES WHICH ARE SHARED BY ALL FRAMES
ctrcoords1=32-ctr:61+ctr 
ctrcoords2=32-ctr:32+ctr 

%%
% >>>>>>>>>>>>>>>>>>>>RECON FROM HERE<<<<<<<<<<<<<<<<<<<<<<<<<<<

% 2: estimate subspaces
nav_parameter_dim1 = squeeze(kspace(ctrcoords1,ctrcoords2,:,:,1));
nav_estimate_1= subspace_estimator_multicoil(nav_parameter_dim1,L3);

nav_parameter_dim2 = squeeze(kspace(ctrcoords1,ctrcoords2,:,1,:));
nav_estimate_2= subspace_estimator_multicoil(nav_parameter_dim2,L4);

% 3: initialize other operators
if strcmp(sparsity_transform,'wavelet')
    Psi=opWavelet2(res,res,'Daubechies') %wavelet operator (uses SPOT toolbox (+ other dependencies maybe?)
elseif strcmp(sparsity_transform,'TV')
    % to do: can be made ~15 times faster with a finite difference operator
%     Psi=opConvolve(res,res,[-1 1],[0 0],'truncated')* opConvolve(res,res,[-1 1]',[0 0],'truncated') %2D TV operator
elseif strcmp(sparsity_transform,'TVcomplex')
%     Psi=TV2op()
    Psi1=opConvolve(res1,res2,[-1 1],[0 0],'truncated')
    Psi2=opConvolve(res1,res2,[-1 1]',[0 0],'truncated')
    Psi=[Psi1;Psi2]
end

sens=squeeze(sens);
sens_normalized=bsxfun(@rdivide,sens,(eps+sqrt(sum(abs(sens).^2,3)))); 
F=MCFopClass;
set_MCFop_Params(F,(sens_normalized),[res1,res2],[tensorsize(4),tensorsize(5)]);

%4 zero-filled recon
P0=F'*kspace;             
P1_0=reshape(P0,unfoldedIsize); %1-unfolding of zero filled recon (what is the shape of this matrix?)
kspace_1=reshape(kspace,unfoldedKsize);

figure(4); imshow(abs(P0(:,:,1,1,1)),[]); axis off; title('zero filled recon of one frame')
figure(5); imshow(angle(P0(:,:,1,1,1)),[]); axis off; title('phase of zero filled recon of one frame')
figure(6); immontage4D(squeeze(abs(P0)));
figure(7); immontage4D(squeeze(angle(P0)),[0 2*pi]);


%% ALGO 


%initialize matrices
[Phi,G,C,A,B,Y,Z]= init_G0(P1_0,Psi,nav_estimate_1,nav_estimate_2,Lg);                    

MSE=[]; 
for iter=1:niter
    MSE=visualize_convergence(iter,MSE,G,C,Phi,[],imagesize,20,33);
    Ak=soft_thresh_A(G,Y,alpha,lambda,Psi);                     %15
    Bk=soft_thresh_B(C,Z,mu,beta);                              %16
    Gk=precon_conj_grad_G(G,C,Ak,Y,alpha,Psi,kspace_1,Phi,F);       %17
    Ck=precon_conj_grad_C(Gk,C,Bk,Z,beta,kspace_1,Phi,F);           %18
    Yk=Y+alpha*(Ak-Psi*Gk);
    Zk=Z+beta.*(Bk-Ck);
    
    G=Gk; C=Ck; Y=Yk; Z=Zk; %update iteration
    
% alpha=alpha*1.5; beta=beta*1.5;    %optional: increasing beta & alpha
end
