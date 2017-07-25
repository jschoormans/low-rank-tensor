clear all; close all; clc
if ispc
%     cd('L:\basic\divi\Ima\parrec\Jasper\Low_Rank_2017_07_11')
    cd('L:\basic\divi\Ima\parrec\Jasper\Low_Rank_2017_07_17')
    addpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations\tensor\low-rank-tensor'))
else
%     cd('/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/Low_Rank_2017_07_11')
%         cd('/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/Low_Rank_2017_07_17')
        cd('/home/qzhang/lood_storage/divi/Ima/parrec/Kerry/Data/2017_07_24_LRT_DTI_T2prep')
        vars_kerry;
    
end

% MR = MRecon('lo_12072017_1926046_17_2_wipvfat2preplowrankV4.raw');
% MR=MRecon('lo_17072017_1603264_2_2_wipvfat2preplowrankV4.raw')
% MR=MRecon('lo_17072017_1635333_4_2_wipvfat2preplowrankV4.raw')
MR=MRecon('lr_24072017_2013158_9_2_wipvfat2preplowrankV4.raw')

%     

MR.Parameter.Recon.RemoveMOversampling='No';
MR.Parameter.Recon.RemovePOversampling='No';
MR.Parameter.Recon.ArrayCompression='Yes';
MR.Parameter.Recon.ACNrVirtualChannels=4;

MR.Parameter.Parameter2Read.typ = 1;
MR.Parameter.Recon.ImmediateAveraging='No'

% MR.Parameter.Recon.CoilCombination = 'pc';
% MR.Parameter.Cardiac.RetroHoleInterpolation = 'No';
% MR_data.Parameter.Cardiac.Synchronization = 'Yes';
% MR.Parameter.Cardiac.RetroPhases = 24;

MR.Parameter.Parameter2Read.Update;

% load data
disp('readdata')
MR.ReadData;
MR.RandomPhaseCorrection;
disp('corrections...')
MR.RemoveOversampling;
MR.PDACorrection;
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;

disp('sortdata')
MR.SortData;
MR.K2IM;
%% 
% K=bart('cc -p3',MR.Data);
data = double(MR.Data);
K=MR.Data(124,:,:,:,:,:,:,:,:,:,:,:);
K = sum(K, 12);
% remove stupid checkerboard pattern
che=create_checkerboard([1,size(K,2),size(K,3)]);
K=bsxfun(@times,K,che);

% a = sum(K(:,:,:,:,1,:,1,1,1,1),6)./sum(K(:,:,:,:,1,:,1,1,1,1)~=0,6);
a = K(:,:,:,:,1);
sens=bart('ecalib -S -m1',a);
du=squeeze(K);

par_dim1 = 9; par_dim2 = 6;
du_sorted = zeros(size(du,1),size(du,2),size(du,3),par_dim1,par_dim2);
for ii = 1:par_dim1
    for jj = 1:par_dim2
        du_sorted(:,:,:,ii,jj) = du(:,:,:,((jj-1) * par_dim1 +ii));
    end
end
mask=squeeze(du(:,:,1,:,:))~=0;

figure(1); clf; immontage4D(mask,[0 1]);
xlabel('Parameter 1'); ylabel('Parameter 2');

% estimate subspaces
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/tensor/low-rank-tensor'))
sparsity_transform='TVcomplex'

res1=size(du,1);
res2=size(du,2);
L3=2;
L4=2;

ctr=3
ctrcoords1=32-ctr:61+ctr %???
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
imshow(abs(Q),[])
clear Q J 


figure(997); 
 Q=[]
for ii=1:size(P_recon,4)
    J=[];
    for jj=1:size(P_recon,5);
        J=[J,abs(P_recon(:,:,1,ii,jj))];
    end
    Q=[Q;J];
end
imshow(abs(Q),[])
clear Q J 