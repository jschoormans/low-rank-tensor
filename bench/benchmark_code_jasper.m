% recon scan 16 (high res)
addpath(genpath('L:\basic\divi\Projects\cosart\tensor\low-rank-tensor'))
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\imagine'))
tic
fprintf('\n == Benchmark == \n',toc)
fprintf('Loading data... \n',toc)

load('benchmark-data.mat')
% load
params.GPU=1;
params.visualization=1; 
params.TVoption=4;
fprintf('Starting recon... \n',toc)
benchstart=tic;
params.GPUdouble=0; %for GPUs supporting double 
params.lambda=2
params.Lg=25;
params.L3=5
params.L4=5

params.mu=15;
params.alpha=20	
params.beta=20
params.niter=20;
params.automu=0
params.increase_penalty_parameters=1
params.TVoption=2
params.G.precon=0

P_recon=LRT_recon_test(kspaceinput,squeeze(sens),params);

imagine(squeeze(P_recon))
fprintf('\n == Benchmark - time LRT_recon : %4.2f == \n\n',toc(benchstart))