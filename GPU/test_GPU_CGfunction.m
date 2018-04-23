% test recon 
cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2018_01_31_VFA');
clear
load('tesdata.mat')
params.GPU=0;
P_recon=LRT_recon(kspaceinput,squeeze(sens),params);


%%
tic
GG=gpuArray(G);
CG=gpuArray(C);
AkG=gpuArray(Ak);
YG=gpuArray(Y);
kspace_1G=gpuArray(kspace_1);
PhiG=gpuArray(Phi);
PsiG=TV_GPU(249,66,2);
toc 

%%
tic;
Gk=precon_conj_grad_G(G,C,Ak,Y,alpha,Psi,kspace_1,Phi,F,params);       %17
tCPU=toc; 

tic;
Gk=precon_conj_grad_G(GG,CG,AkG,YG,alpha,PsiG,kspace_1G,PhiG,F,params);       %17
tGPU=toc; 

fprintf('CG G: time taken CPU: %4.2f seconds \n', tCPU)
fprintf('CG G: time taken GPU: %4.2f seconds \n', tGPU)
fprintf('speed up: %4.2f times \n', tCPU/tGPU)


