function P_recon=LRT_recon(kspace,sens,params)
% 4D Low-Rank Tensor reconstruction 

% recon is run as: P_recon=LRT_recon(kspace,sens,params)
% params is a struct- should be made with params_init

% inputs:
% kspace: undersampled kspace with dimensions (x,y,coils,param1,param2)
% sens  : sens maps - size (x,y,coils) 
% params: params is a struct- should be made with params_init

% based on IEEE-TMI paper by Jingfei He  ,2016 

% 2017 J Schoormans and Q Zhang - AMC Amsterdam

fprintf('------------------------------ \n')
fprintf('Low-Rank Tensor Reconstruction \n')
fprintf('------------------------------ \n')
fprintf('   J Schoormans & Q Zhang  \n')
fprintf('       Amsterdam 2017      \n \n')



%%
assert(params.Lg<=params.L3*params.L4,'reducte spatial rank!'); 

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


% >>>>>>>>>>>>>>>>>>>>RECON FROM HERE<<<<<<<<<<<<<<<<<<<<<<<<<<<
[Kx1,Ky1,Kx2,Ky2]=findSharedKpoints(mask);

% 2: estimate subspaces: generalized for non-square shared k-points
nav_parameter_dim1=[];
for iter=1:length(Kx1)
    nav_parameter_dim1=cat(1,nav_parameter_dim1,(kspace(Kx1(iter),Ky1(iter),:,:,1)));
end
[nav_estimate_1,eigenvals_1]= subspace_estimator_multicoil(squeeze(nav_parameter_dim1),params.L3);

nav_parameter_dim2=[];
for iter=1:length(Kx2)
    nav_parameter_dim2=cat(1,nav_parameter_dim2,(kspace(Kx2(iter),Ky2(iter),:,1,:)));
end
[nav_estimate_2,eigenvals_2]= subspace_estimator_multicoil(squeeze(nav_parameter_dim2),params.L4);

% 3: initialize other operators
if strcmp(params.sparsity_transform,'wavelet')
    Psi=opWavelet2(res1,res2,'Daubechies'); %wavelet operator (uses SPOT toolbox (+ other dependencies maybe?)
    operatorsize=[96,96]; % not sure - to do
elseif strcmp(params.sparsity_transform,'TV')
    Psi1=opConvolve(res1,res2,[-1 1],[0 0],'cyclic');
    Psi2=opConvolve(res1,res2,[-1 1]',[0 0],'cyclic');
    Psi=[Psi1;Psi2];
    operatorsize=[res1,res2*2];
elseif strcmp(params.sparsity_transform,'TVOP')
    Psi=TVOP;
    operatorsize=[res1,res2*2]; % not sure - to do 
elseif strcmp(params.sparsity_transform,'I')
    Psi=opEye(res1*res2);
    operatorsize=[res1,res2]; % not sure - to do 

else
    error('sparsity_transform not recognized')
end

sens_normalized=bsxfun(@rdivide,sens,(eps+sqrt(sum(abs(sens).^2,3)))); 
F=MCFopClass;
set_MCFop_Params(F,(sens_normalized),[res1,res2],[tensorsize(4),tensorsize(5)]);

%4 zero-filled recon
P0=F'*kspace;             
P1_0=reshape(P0,unfoldedIsize); %1-unfolding of zero filled recon (what is the shape of this matrix?)

if params.scaleksp
    [kspace,scaling]= scaleksp(kspace,P0); % scale kspace to ensure consistency over params;
    params.Imref=params.Imref./scaling; %scale ref image with same scaling;
end 

kspace_1=reshape(kspace,unfoldedKsize);

figure(4); imshow(abs(P0(:,:,1,1,1)),[]); axis off; title('zero filled recon of one frame')
figure(5); imshow(angle(P0(:,:,1,1,1)),[]); axis off; title('phase of zero filled recon of one frame')
figure(6); immontage4D(squeeze(abs(P0)));
figure(7); immontage4D(squeeze(angle(P0)),[-pi pi]);
set(0,'DefaultAxesColorOrder',jet(max([size(nav_estimate_1,2), size(nav_estimate_2,2)]))); 
figure(8); plot(abs(nav_estimate_1)); colorbar
figure(9); plot(abs(nav_estimate_2)); colorbar
figure(10); hold on; plot(eigenvals_1./max(eigenvals_1(:)),'r'); plot(params.L3,eigenvals_1(params.L3)./max(eigenvals_1(:)),'ro');...
    plot(eigenvals_2./max(eigenvals_2(:)),'b');plot(params.L4,eigenvals_2(params.L4)./max(eigenvals_2(:)),'bo') ;hold off;
title('eigenvalues for the two subspaces (1=red,2-blue)')
%% ALGO 
alpha=params.alpha;
beta=params.beta; 
mu=params.mu; 
lambda=params.lambda; 

%initialize matrices
if params.inspectLg
params.Lg=params.L3*params.L4;end

[Phi,G,C,A,B,Y,Z]= init_G0(P1_0,Psi,nav_estimate_1,nav_estimate_2,params.Lg);  

if params.inspectLg;
    C_energies=sum(abs(C).^2,2);
    C_energies=C_energies./max(C_energies(:));
    figure(11); plot(C_energies,'ko-');title('rank relative energies for C');
    fprintf('relE %f: \n',C_energies)
    params.Lg=input('Choose spatial rank: ');
    [Phi,G,C,A,B,Y,Z]= init_G0(P1_0,Psi,nav_estimate_1,nav_estimate_2,params.Lg);  

end


MSE=[]; 

sensmask=sum(sens_normalized,3)~=0;
sensmask=reshape(sensmask,[res1*res2 1]);

for iter=1:params.niter
    fprintf('\n Outer iteration %i of %i \n',iter,params.niter)
    MSE=visualize_convergence(iter,MSE,G,C,Phi,params.Imref,imagesize,params.x,params.y);
    Ak=soft_thresh_A(G,Y,alpha,lambda,Psi,operatorsize);                     %15
    Bk=soft_thresh_B(C,Z,mu,beta);                              %16
    Gk=precon_conj_grad_G(G,C,Ak,Y,alpha,Psi,kspace_1,Phi,F,params);       %17
    Ck=precon_conj_grad_C(Gk,C,Bk,Z,beta,kspace_1,Phi,F,params);           %18
    Yk=Y+alpha*(Ak-Psi*Gk);
    Zk=Z+beta.*(Bk-Ck);
    
    G=Gk; C=Ck; Y=Yk; Z=Zk; %update iteration

    if params.increase_penalty_parameters
    alpha=alpha*1.5; beta=beta*1.5; end;   
end

P_recon=G*C*Phi;
P_recon=reshape(P_recon,imagesize);

