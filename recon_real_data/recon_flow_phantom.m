clear all; 
close all; 
clc
% cd('/home/jschoormans/lood_storage/divi/Temp/Lukas/reconstruction/test-data/spiralshape')

if ispc()
cd('L:\basic\divi\Temp\Lukas\reconstruction\test-data\spiralshape')
else
    cd('/home/jschoormans/lood_storage/divi/Temp/Lukas/reconstruction/test-data/spiralshape')
end
file='cs_12012017_1913584_9_2_wipk64r6mode1np50t3flip1senseV4.raw'
% file='cs_12012017_1905246_8_2_wipk64r6mode1np50t3flip0senseV4.raw'
%% copied from lukas' code'

mrecon=MRecon(file)
mrecon.Parameter.Recon.CoilCombination = 'pc';

% cardiac binning
mrecon.Parameter.Cardiac.RetroHoleInterpolation = 'No';
mrecon.Parameter.Cardiac.Synchronization = 'Retrospective';
mrecon.Parameter.Cardiac.RetroPhases = 24;

mrecon.Parameter.Parameter2Read.typ = 1;
mrecon.Parameter.Parameter2Read.Update;

% load data
mrecon.ReadData;

mrecon.RandomPhaseCorrection;
mrecon.RemoveOversampling;
mrecon.PDACorrection;
mrecon.DcOffsetCorrection;
mrecon.MeasPhaseCorrection;
mrecon.SortData;

%%
mrecon.K2IM
%%
K=mrecon.Data(32,:,:,:,:,:,:,:,:,:,:,:);
size(K)
% remove stupid checkerboard pattern
che=create_checkerboard([1,64,64]);
K=bsxfun(@times,K,che);

a = sum(K(:,:,:,:,1,:,1,1,1,1),6)./sum(K(:,:,:,:,1,:,1,1,1,1)~=0,6);
sens=bart('ecalib -r 20 -m1',a);

%
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
sparsity_transform='TVcomplex'

res1=size(du,1);
res2=size(du,2);
L3=2;
L4=2;

ctr=1
ctrcoords1=32-ctr:32+ctr %???
ctrcoords2=32-ctr:32+ctr
% >>>>>>>>>>>>>>>>>>>>RECON FROM HERE<<<<<<<<<<<<<<<<<<<<<<<<<<<

% 2: estimate subspaces
nav_parameter_dim1 = squeeze(du(ctrcoords1,ctrcoords2,:,:,1));
nav_estimate_1= subspace_estimator_multicoil(nav_parameter_dim1,L3);

nav_parameter_dim2 = squeeze(du(ctrcoords1,ctrcoords2,:,1,:));
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
    Psi1=opConvolve(res1,res2,[-1 1],[0 0],'truncated')
    Psi2=opConvolve(res1,res2,[-1 1]',[0 0],'truncated')
    Psi=[Psi1;Psi2]
end

sens=squeeze(sens);
sens_normalized=bsxfun(@rdivide,sens,(eps+sqrt(sum(abs(sens).^2,3)))); 
F=MCFopClass;
set_MCFop_Params(F,(sens_normalized),[res1,res2],[tensorsize(4),tensorsize(5)]);

%4 zero-filled recon
P0=F'*du;             
P1_0=reshape(P0,unfoldedIsize); %1-unfolding of zero filled recon (what is the shape of this matrix?)
du_1=reshape(du,unfoldedKsize);

figure(4); imshow(abs(P0(:,:,1,1,1)),[]); axis off; title('zero filled recon of one frame')
figure(5); imshow(angle(P0(:,:,1,1,1)),[]); axis off; title('phase of zero filled recon of one frame')
figure(6); immontage4D(squeeze(abs(P0)));
figure(7); immontage4D(squeeze(angle(P0)),[0 2*pi]);

%% 
% lowres_phase_estimate=exp(-1i*angle(P0(:,:,1,1,1)));
% Psi=Psi*opDiag(lowres_phase_estimate(:))

%% ALGO 
%initialize parameters
alpha= 2;         %penalty parameter >0
beta=  2;         %penalty parameter >0
alpha=0.01/(mean(abs(du(du~=0))));
beta=alpha; %from paper 

lambda=1e-8;        %sparsity parameter
mu=1e-1;           %sparsity parameter

% Lg=L3*L4;             %rank of spatial dimension
Lg=3;
niter=10;

%initialize matrices
[Phi,G,C,A,B,Y,Z]= init_G0(P1_0,Psi,nav_estimate_1,nav_estimate_2,Lg);                    

MSE=[]; 
for iter=1:niter
    MSE=visualize_convergence(iter,MSE,G,C,Phi,[],imagesize,20,33);
    Ak=soft_thresh_A(G,Y,alpha,lambda,Psi);                     %15
    Bk=soft_thresh_B(C,Z,mu,beta);                              %16
    Gk=precon_conj_grad_G(G,C,Ak,Y,alpha,Psi,du_1,Phi,F);       %17
    Ck=precon_conj_grad_C(Gk,C,Bk,Z,beta,du_1,Phi,F);           %18
    Yk=Y+alpha*(Ak-Psi*Gk);
    Zk=Z+beta.*(Bk-Ck);
    
    G=Gk; C=Ck; Y=Yk; Z=Zk; %update iteration
    
% alpha=alpha*1.5; beta=beta*1.5;    %optional: increasing beta & alpha
end