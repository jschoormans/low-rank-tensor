function [P_recon,params]=LRT_recon_test(kspace,sens,params)
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
assert(params.Lg<=params.L3*params.L4,'reduce spatial rank!'); 
assert(~xor(isempty(params.nav_estimate_1),isempty(params.nav_estimate_2)),'Input either both or no subspaces!')
res1=size(kspace,1);
res2=size(kspace,2);
tensorsize=size(kspace);
imagesize=tensorsize; imagesize(3)=1; 
unfoldedIsize=[size(kspace,1)*size(kspace,2),size(kspace,4)*size(kspace,5)];                %coil combined
unfoldedKsize=[size(kspace,1)*size(kspace,2)*size(kspace,3),size(kspace,4)*size(kspace,5)];     %coils separate

%%

if params.GPU; 
    kspace=gpuArray(kspace); 
end

mask=squeeze(kspace(:,:,1,:,:))~=0;

% >>>>>>>>>>>>>>>>>>>>RECON FROM HERE<<<<<<<<<<<<<<<<<<<<<<<<<<<

if isempty(params.nav_estimate_1)                  % subspace estimation 
    [Kx1,Ky1,Kx2,Ky2]=findSharedKpoints(mask,params);
    fprintf('Estimating subspaces...\n')
    fprintf('# of shared ky-kz points of dim 1: %i \n',numel(Kx1))
    fprintf('# of shared ky-kz points of dim 2: %i \n',numel(Kx2))
    
    % 2: estimate subspaces: generalized for non-square shared k-points
    nav_parameter_dim1=[];
    for iter=1:length(Kx1)
        nav_parameter_dim1=cat(1,nav_parameter_dim1,(kspace(Kx1(iter),Ky1(iter),:,:,params.subspacedim1)));
    end
    [nav_estimate_1,params.eigenvals_1]= subspace_estimator_multicoil(squeeze(nav_parameter_dim1),params.L3);
    
    nav_parameter_dim2=[];
    for iter=1:length(Kx2)
        nav_parameter_dim2=cat(1,nav_parameter_dim2,(kspace(Kx2(iter),Ky2(iter),:,params.subspacedim2,:)));
    end
    [nav_estimate_2,params.eigenvals_2]= subspace_estimator_multicoil(squeeze(nav_parameter_dim2),params.L4);
    
    params.nav_estimate_1=nav_estimate_1;
    params.nav_estimate_2=nav_estimate_2;
    params.params.eigenvals_1=params.eigenvals_1;
    params.params.eigenvals_2=params.eigenvals_2;
else
        fprintf('Subspaces user-defined...\n')
end

% 3: initialize other operators
if strcmp(params.sparsity_transform,'wavelet')
    Psi=opWavelet2(res1,res2,'Daubechies'); %wavelet operator (uses SPOT toolbox (+ other dependencies maybe?)
    operatorsize=[96,96]; % not sure - to do
elseif strcmp(params.sparsity_transform,'TV')
    Psi=TV_GPU(res1,res2,params.TVoption,params.GPUdouble);
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

if params.normalize_sense %find out how it should be done...
    sensvec=reshape(sens,[size(sens,1)*size(sens,2),size(sens,3)]);
    sens1_mag = reshape(vecnorm(reshape(sensvec, [], size(sens,3)).'), [size(sens,1),size(sens,2)]);
    sens_normalized = bsxfun(@rdivide, sens, sens1_mag);
    sens_normalized(find(isnan(sens_normalized))) = 0;
else sens_normalized=sens;
end
F=MCFopClass;
if params.GPU
    if params.GPUdouble
        set_MCFop_Params(F,gpuArray(double(sens_normalized)),[res1,res2],[tensorsize(4),tensorsize(5)]);
    else
        set_MCFop_Params(F,gpuArray(single(sens_normalized)),[res1,res2],[tensorsize(4),tensorsize(5)]);
    end
else
    set_MCFop_Params(F,(sens_normalized),[res1,res2],[tensorsize(4),tensorsize(5)]);
end

%4 zero-filled recon


P0=F'*(kspace);
P1_0=reshape(P0,unfoldedIsize); %1-unfolding of zero filled recon (what is the shape of this matrix?)

if params.scaleksp
    [kspace,scaling]= scaleksp(kspace,P0); % scale kspace to ensure consistency over params;
    params.Imref=params.Imref./scaling; %scale ref image with same scaling;
    P0=F'*(kspace);             
    P1_0=reshape(P0,unfoldedIsize); %1-unfolding of zero filled recon (what is the shape of this matrix?)
end 

kspace_1=reshape(kspace,unfoldedKsize);

if params.visualization==1
figure(21); subplot(211);   imshow(abs(P0(:,:,1,params.subspacedim2,params.subspacedim1)),[]); axis off; title('zero filled recon of one frame')
figure(21); subplot(212);   imshow(angle(P0(:,:,1,params.subspacedim2,params.subspacedim1)),[]); axis off; title('phase of zero filled recon of one frame')
figure(22);subplot(311);    immontage4D(mask,[0 1]); xlabel('Parameter 1'); ylabel('Parameter 2');
figure(22); subplot(312);   immontage4D(squeeze(abs(P0)));
figure(22); subplot(313);   immontage4D(squeeze(angle(P0)),[-pi pi]);
set(0,'DefaultAxesColorOrder',jet(max([size(params.nav_estimate_1,2), size(params.nav_estimate_2,2)]))); 
figure(23); subplot(221); cla;plot(abs(params.nav_estimate_1)); colorbar
figure(23); subplot(222); cla;plot(abs(params.nav_estimate_2)); colorbar
figure(23); subplot(212); cla;hold on; plot(params.eigenvals_1./max(params.eigenvals_1(:)),'r'); plot(params.L3,params.eigenvals_1(params.L3)./max(params.eigenvals_1(:)),'ro');...
    plot(params.eigenvals_2./max(params.eigenvals_2(:)),'b');plot(params.L4,params.eigenvals_2(params.L4)./max(params.eigenvals_2(:)),'bo') ;hold off;
title('eigenvalues for the two subspaces (1=red,2-blue)');
drawnow;
end;
%% ALGO 
alpha=params.alpha;
beta=params.beta; 
mu=params.mu; 
lambda=params.lambda; 

%initialize matrices
if params.inspectLg %give temporary spatial rank (max rank)
params.Lg=params.L3*params.L4;end

[Phi,G,C,Ak,Bk,Y,Z]= init_G0(P1_0,Psi,params.nav_estimate_1,params.nav_estimate_2,params.Lg);  

if params.inspectLg; % option to check energies of spatial ranks before choosing rank
    C_energies=sum(abs(C).^2,2);
    C_energies=C_energies./max(C_energies(:));
    if params.visualization;    figure(11); plot(C_energies,'ko-');title('rank relative energies for C'); end
    fprintf('relE %f: \n',C_energies)
    params.Lg=input('Choose spatial rank: ');
    [Phi,G,C,Ak,Bk,Y,Z]= init_G0(P1_0,Psi,params.nav_estimate_1,params.nav_estimate_2,params.Lg);  
end

if params.GPU
    % convert all to singles (test)
    if ~params.GPUdouble
        G=single(G);
        C=single(C);
        Ak=single(Ak);
        Bk=single(Bk);
        Y=single(Y);
        Z=single(Z);
        kspace_1=single(kspace_1);
        Phi=single(Phi);
    else
        G=double(G);
        C=double(C);
        Ak=double(Ak);
        Bk=double(Bk);
        Y=double(Y);
        Z=double(Z);
        kspace_1=double(kspace_1);
        Phi=double(Phi);
    end
    % define GPU Arrays 
    G=gpuArray(G);
    C=gpuArray(C);
    Ak=gpuArray(Ak);
    Bk=gpuArray(Bk);
    Y=gpuArray(Y);
    Z=gpuArray(Z);
    kspace_1=gpuArray(kspace_1);
    Phi=gpuArray(Phi);
end

MSE=[]; 
for iter=1:params.niter
    params.iter=iter; 
    fprintf('\n=== Outer iteration %i of %i === \n',params.iter,params.niter)
    fprintf('alpha = %4.2f | beta = %4.2f \n',alpha,beta)

    if params.visualization;
    MSE=visualize_convergence(params.iter,MSE,G,C,Phi,params.Imref,imagesize,params.x,params.y,kspace_1,F,Psi,Ak); end
    
    Aold=Ak; Bold=Bk;
    [Ak,lambda]=soft_thresh_A(G,Y,alpha,lambda,Psi,operatorsize,params);                     %15
    [Bk,mu]=soft_thresh_B(C,Z,mu,beta,params);                              %16
    
    %new addition:
    Y=Y+alpha*(Ak-(Psi*G));
    Z=Z+beta.*(Bk-C);  
    
    Gk=conj_grad_G_new(G,C,Ak,Y,alpha,Psi,kspace_1,Phi,F,params);       %17
    Ck=conj_grad_C_new(Gk,C,Bk,Z,beta,kspace_1,Phi,F,params);           %18
    Yk=Y+alpha*(Ak-(Psi*Gk));
    Zk=Z+beta.*(Bk-Ck);
    
    G=Gk; C=Ck; Y=Yk; Z=Zk; %update iteration
    
    % temp norms
    rnorm1(iter)=norm(Ak-(Psi*G));
    rnorm2(iter)=norm(Bk-C);
    
    snorm1(iter)=norm(-(alpha/2)*(Ak-Aold));
    snorm2(iter)=norm(-(beta/2)*(Bk-Bold));

    
    [obj(iter),l2(iter),l1(iter),l12(iter)]=objective(kspace_1,F,G,C,Phi,Ak,Bk,lambda,mu);
    
    figure(2000);
    subplot(311); hold on; 
    semilogy(obj,'k.-');
    semilogy(l1,'r.');
    semilogy(l2,'g.');
    semilogy(l12,'b.'); hold off; 
    title('objective'); xlabel('iterations');
    legend('total objective','l1','l2','l12')
    subplot(323); semilogy(rnorm1,'k.-'); title('r-norm 1'); xlabel('iterations'); 
    subplot(324); semilogy(rnorm2,'k.-'); title('r-norm 2'); xlabel('iterations'); drawnow; 
    subplot(325); semilogy(snorm1,'k.-'); title('s-norm 1'); xlabel('iterations'); 
    subplot(326); semilogy(snorm2,'k.-'); title('s-norm 2'); xlabel('iterations'); drawnow; 

    
    
    if params.increase_penalty_parameters
    alpha=alpha*1.5; beta=beta*1.5; end;   
end

% after last iteration
if params.visualization;
    MSE=visualize_convergence(params.iter,MSE,G,C,Phi,params.Imref,imagesize,params.x,params.y,kspace_1,F,Psi,Ak); end

P_recon=G*C*Phi;
P_recon=reshape(P_recon,imagesize);
end

function [obj,l2,l1,l12]=objective(d,F,G,C,Phi,A,B,lambda,mu)
l2=norm(d-F*(G*(C*Phi)));
l1=mu*l1norm(B);
l12=lambda*l12norm(A);
obj=l2+l1+l12;
end

function l=l1norm(B)
l=sum(abs(B(:)));
end

function l=l12norm(A)
l=sum((sum(abs(A).^2,2)).^(1/2)); % col or row??
end




