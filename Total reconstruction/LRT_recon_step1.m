function [P_recon,goldstandardscaled,C,G1,G2,G3] = LRT_recon_step1(kspace,sens)

% Define parameters
params=params_init(kspace);
%

[P_recon,goldstandardscaled,C,G1,G2,G3]=LRT_recon_step2(kspace,sens,params);
% visualize recon %IS DIFFERENT THAN IN RECON??
% figure(1000); immontage4D(squeeze(abs(P_recon)),[]); % line 30 to 40 once indented on 15-11-2017
% figure(1001); immontage4D(squeeze(angle(P_recon)),[-pi pi]);
% figure(1002); imshow(angle(P_recon(:,:,1,1,1)),[-pi pi]);
% %
% cplkdiff=squeeze((P_recon(:,:,:,:,1))./(P_recon(:,:,:,:,4)));
% 
%   figure(1003); immontage4D(reshape(angle(cplkdiff),[192 192 6 4]),[-pi/8 pi/8])
% %   figure(1003); immontage4D(reshape(angle(cplkdiff),[160 160 6 5]),[-pi/8 pi/8])
% 
% % cplkdiff=squeeze((P_recon(:,:,:,:,1))./(P_recon(:,:,:,:,4)));
% %   figure(1003); immontage4D(reshape(angle(cplkdiff),[160 160 6 5]),[-pi/8 pi/8])



