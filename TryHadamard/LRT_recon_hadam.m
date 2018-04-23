function [P_recon,goldstandardscaled,C,G,nav_estimate_1,nav_estimate_2]=LRT_recon_hadam(kspace,sens,params)
% 4D Low-Rank Tensor reconstruction 

% recon is run as: P_recon=LRT_recon_hadam(kspace,sens)

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
close all
sens=sens+1e-4;
sens=squeeze(sens);

assert(~xor(isempty(params.nav_estimate_1),isempty(params.nav_estimate_2)),'Input either both or no subspaces!')
res1=size(kspace,1);
res2=size(kspace,2);
tensorsize=size(kspace);
imagesize=tensorsize; imagesize(3)=1; 
unfoldedIsize=[size(kspace,1)*size(kspace,2),size(kspace,4)*size(kspace,5)];                %coil combined
unfoldedKsize=[size(kspace,1)*size(kspace,2)*size(kspace,3),size(kspace,4)*size(kspace,5)];     %coils separate

%%
% FIND MASK
mask=squeeze(kspace(:,:,1,:,:))~=0;

if params.hadamard == 1;
% do hadamard of kspace

% we have to use kspace_nav further on because we need a different
    % quantity below. Todo: use matrixform as operator.
    kspace_nav(:,:,:,:,1)=kspace(:,:,:,:,1)+kspace(:,:,:,:,2)+kspace(:,:,:,:,3)+kspace(:,:,:,:,4);
    kspace_nav(:,:,:,:,2)=kspace(:,:,:,:,1)+kspace(:,:,:,:,2)-kspace(:,:,:,:,3)-kspace(:,:,:,:,4);
    kspace_nav(:,:,:,:,3)=kspace(:,:,:,:,1)-kspace(:,:,:,:,2)+kspace(:,:,:,:,3)-kspace(:,:,:,:,4);
    kspace_nav(:,:,:,:,4)=kspace(:,:,:,:,1)-kspace(:,:,:,:,2)-kspace(:,:,:,:,3)+kspace(:,:,:,:,4);
    kspace_nav=0.5*kspace_nav;
else kspace_nav = kspace;
end

% >>>>>>>>>>>>>>>>>>>>RECON FROM HERE<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Changed by Bobby on 18-12-2017 to add columns and rows in navigation.
if isempty(params.nav_estimate_1)                  % subspace estimation 
    nav_parameter_dim1=[];
    nav_parameter_dim2=[];
    for column = 1:size(params.columns,2);
        [Kx1,Ky1]=findSharedKpoints_param1(mask,params.columns(column));
        fprintf('Estimating subspaces...\n')
        fprintf('# of shared ky-kz points of dim 1: %i \n',numel(Kx1))
        
        % 2: estimate subspaces (1): generalized for non-square shared k-points
        for iter=1:length(Kx1)
            nav_parameter_dim1=cat(1,nav_parameter_dim1,(kspace_nav(Kx1(iter),Ky1(iter),:,:,params.columns(column))));
    %         if iter == 2000;
    %             break
    %         end
        end
    end
    [nav_estimate_1,params.eigenvals_1]= subspace_estimator_multicoil(squeeze(nav_parameter_dim1),params.L3);
    
    %add rows here
    for row = 1:size(params.rows,2);
        [Kx2,Ky2]=findSharedKpoints_param2(mask,params.rows(row));
        fprintf('# of shared ky-kz points of dim 2: %i \n',numel(Kx2))
        
         % 2: estimate subspaces (2): generalized for non-square shared k-points
        for iter=1:length(Kx2)
            nav_parameter_dim2=cat(1,nav_parameter_dim2,(kspace_nav(Kx2(iter),Ky2(iter),:,params.rows(row),:)));
    %         if iter == 2000;
    %             break
    %         end
        end
    end
    [nav_estimate_2,params.eigenvals_2]= subspace_estimator_multicoil2(squeeze(nav_parameter_dim2),params.L4);
    
    %Added by Bobby 16-01-2018 to test with subspace from 'tucker'
%     nav_estimate_1 = [];
%     nav_estimate_2 = [];
%     facmats = open('/home/barunderkamp/lood_storage/divi/Projects/4dflow_lrt/data/2017_11_28_FlowPhantomFully/recon_out/tucker_ownlrtcompare_09-01-2018/factormatricesfullspatialrankwithoversampling.mat');
%     nav_estimate_1 = facmats.facmats{3};
%     nav_estimate_2 = facmats.facmats{4};
    
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

if params.normalize_sense %find out how it should be done...
    sensvec=reshape(sens,[size(sens,1)*size(sens,2),size(sens,3)]);
    sens1_mag = reshape(vecnorm(reshape(sensvec, [], size(sens,3)).'), [size(sens,1),size(sens,2)]);
    sens_normalized = bsxfun(@rdivide, sens, sens1_mag);
    sens_normalized(find(isnan(sens_normalized))) = 0;
else sens_normalized=sens;
end
F=MCFopClass;
% Modified by Bobby 16-01-2018
% set_MCFop_Params(F,(sens_normalized),[res1,res2],[tensorsize(4),tensorsize(5)]);
 set_MCFop_Params(F,(sens_normalized),[res1,res2],[size(kspace,4),size(kspace,5)]);
%4 zero-filled recon


P0=F'*(kspace);   

if params.scaleksp
    [kspace,scaling]= scaleksp(kspace,P0); % scale kspace to ensure consistency over params;
    params.Imref=params.Imref./scaling; %scale ref image with same scaling;
    goldstandardscaled=params.Imref;
    P0=F'*(kspace);             
end 
fprintf('test nu even hoe de scaling is en of 0.1 goeie keus is. \n')
if params.nullbackground;
P0 = (abs(P0)>0.1).*P0;
kspace = F*(P0);
end

if params.hadamard == 1; % we have to use P0_2 further on because we need a different
    % quantity below. Todo: use matrixform as operator.
    P0_2(:,:,:,:,1)=P0(:,:,:,:,1)+P0(:,:,:,:,2)+P0(:,:,:,:,3)+P0(:,:,:,:,4);
    P0_2(:,:,:,:,2)=P0(:,:,:,:,1)+P0(:,:,:,:,2)-P0(:,:,:,:,3)-P0(:,:,:,:,4);
    P0_2(:,:,:,:,3)=P0(:,:,:,:,1)-P0(:,:,:,:,2)+P0(:,:,:,:,3)-P0(:,:,:,:,4);
    P0_2(:,:,:,:,4)=P0(:,:,:,:,1)-P0(:,:,:,:,2)-P0(:,:,:,:,3)+P0(:,:,:,:,4);
    P0_2=0.5*P0_2;
else P0_2 = P0;
end


P1_0=reshape(P0_2,unfoldedIsize); %1-unfolding of zero filled recon (what is the shape of this matrix?)

kspace_1=reshape(kspace,unfoldedKsize);

if params.visualize;
% line 118 to 128 once indented on 15-11-2017, unindented 30-11-2017
    figure(21); subplot(211); imshow(abs(P0_2(:,:,1,8,2)),[]); axis off; title('zero filled recon of one frame')
    figure(21); subplot(212);  imshow(angle(P0_2(:,:,1,8,2)),[]); axis off; title('phase of zero filled recon of one frame')
    figure(22);subplot(311); immontage4D(mask,[0 1]); xlabel('Parameter 1'); ylabel('Parameter 2');
    figure(22); subplot(312);  immontage4D(squeeze(abs(P0_2)));
    figure(22); subplot(313); immontage4D(squeeze(angle(P0_2)),[-pi pi]);
    set(0,'DefaultAxesColorOrder',jet(max([size(params.nav_estimate_1,2), size(params.nav_estimate_2,2)]))); 
    figure(23); subplot(221); plot(abs(params.nav_estimate_1)); colorbar
    figure(23); subplot(222); plot(abs(params.nav_estimate_2)); colorbar
    figure(23); subplot(212); hold on; plot(params.eigenvals_1./max(params.eigenvals_1(:)),'r'); plot(params.L3,params.eigenvals_1(params.L3)./max(params.eigenvals_1(:)),'ro');...
        plot(params.eigenvals_2./max(params.eigenvals_2(:)),'b');plot(params.L4,params.eigenvals_2(params.L4)./max(params.eigenvals_2(:)),'bo') ;hold off;
    title('eigenvalues for the two subspaces (1=red,2-blue)');
    drawnow;
end
%% ALGO 
alpha=params.alpha;
beta=params.beta; 
mu=params.mu; 
lambda=params.lambda; 

%initialize matrices
if params.inspectLg
params.Lg=params.L3*params.L4;end

[Phi,G,C,A,B,Y,Z]= init_G0(P1_0,Psi,params.nav_estimate_1,params.nav_estimate_2,params.Lg);  

if params.inspectLg;
    C_energies=sum(abs(C).^2,2);
    C_energies=C_energies./max(C_energies(:));
    figure(11); plot(C_energies,'ko-');title('rank relative energies for C');
    fprintf('relE %f: \n',C_energies)
    params.Lg=input('Choose spatial rank: ');
    [Phi,G,C,A,B,Y,Z]= init_G0(P1_0,Psi,params.nav_estimate_1,params.nav_estimate_2,params.Lg);  
end


MSE=[]; 
for iter=1:params.niter
    params.iter=iter; 
    fprintf('\n Outer iteration %i of %i \n',params.iter,params.niter)
    if params.visualize == 1;
        MSE=visualize_convergence(params.iter,MSE,G,C,Phi,params.Imref,imagesize,params.x,params.y);
    end
    [Ak,lambda]=soft_thresh_A(G,Y,alpha,lambda,Psi,operatorsize,params);                     %15
    [Bk,mu]=soft_thresh_B(C,Z,mu,beta,params);                              %16
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

P_recon = squeeze(P_recon);

if params.hadamard == 1;
    % Prepare transform from hadamard to image domain and resort back.
    psi_hadam = 0.5*[[1,1,1,1];[1,1,-1,-1];[1,-1,1,-1];[1,-1,-1,1]];
    reshape_p_1 = [size(P_recon,1)*size(P_recon,2)*size(P_recon,3),size(P_recon,4)];
    P_recon = reshape(((psi_hadam*((reshape(P_recon,reshape_p_1)).')).'),size(P_recon));
    %
end

