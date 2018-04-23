function P_recon = LRT_recon_1(kspace,sens)

sens=sens+1e-4;
close all;
params=params_init(kspace);

params.Lg=6;
params.inspectLg=false

venc_directions = 4;

params.L3=6; % ranknumber in dimension 3 (cardiac)
params.L4=venc_directions; % rank number in dim 4 (v_enc)

% Added/modified by Bobby 18-12-2017
% params.subspacedim1=2 %kolom for which you calculate SVD 
% params.subspacedim2=8 %row 
params.columns=1:venc_directions; % max number of columns
params.rows=1:24;  % max number of rows
% params.columns=2; % max number of columns
% params.rows=8;  % max number of rows
%

params.sparsity_transform='TV';
params.Imref=[];
params.x=35;
params.y=50;
% params.x=80;
% params.y=16;
% params.y = 80;
params.lambda=8e-1
params.mu=2e0
% sens(sens==0)=1e-2;

params.increase_penalty_parameters=false
params.G.precon=false;
params.G.maxiter=8
params.niter= 20 %initially 20


P_recon=LRT_recon(kspace,squeeze(sens),params);
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



