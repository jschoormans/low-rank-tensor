% recon scan 16 (high res)
addpath(genpath('L:\basic\divi\Projects\cosart\tensor\low-rank-tensor'))
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
params.lambda=0.2
params.Lg=5;
params.mu=5;
params.alpha=2	
params.beta=2
params.niter=5;
params.automu=0
params.increase_penalty_parameters=0
params.TVoption=2
params.G.precon=0

P_recon=LRT_recon_test(kspaceinput,squeeze(sens),params);

fprintf('\n == Benchmark - time LRT_recon : %4.2f == \n\n',toc(benchstart))