% recon simulation
% draft reconstruction code 
clear all; close all; clc;

% 1: make data (most settings in other .m file for now)
uf=0.02; % undersampling factor (excluding center)
noiselevel=0;
ncoils=2;
complexsim=0;
center_for_all_frames=1;

simulation_typ = 'DTI_n_T2W'; %choose from: DWI_n_T2W; VFA_n_T2W; DTI_n_T2W
run create_undersampled_measurement_S.m    
%% 

params=params_init();
params.Lg=1;
params.L3=4; 
params.L4=4;
params.mu=1e4
params.sparsity_transform='TV'
params.Imref=I;
params.x=50
params.y=60;

LRT_recon(du,squeeze(sens),params)