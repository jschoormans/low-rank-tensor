% test recon 
cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2018_01_31_VFA');
clear
load('tesdata.mat')
params.lambda=1e-1;
params.GPU=1; 
params.subspacedim1=2;
params.G.maxiter=40;
params.Lg=12;
params.L3=4;
params.niter=3;
P_recon=LRT_recon(kspaceinput,squeeze(sens),params);
%%
figure(1000); immontage4D(squeeze(abs(P_recon)),[]);

