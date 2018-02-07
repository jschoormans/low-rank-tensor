% recon scan 16 (high res) 
clear; 
cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2018_02_01_VFA')
load('testdata_scan16_sl-200.mat')
% load
params.G.maxiter=10;
params.L3=4;
params.L4=4;
params.subspacedim1=2

params.Lg=6
params.lambda=16e-2;
params.niter=4
params.alpha=2
params.GPUoptions.TVoption=4

P_recon=LRT_recon(kspaceinput(:,:,:,:,:),squeeze(sens(:,:,:,:)),params);
