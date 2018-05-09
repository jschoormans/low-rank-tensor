% recon scan 16 (high res)
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
params.lambda=5e-4
params.Lg=3;
params.mu=0.9;
params.alpha=1	
params.beta=1
params.niter=25;
params.automu=0
params.increase_penalty_parameters=0
P_recon=LRT_recon(kspaceinput,squeeze(sens),params);

fprintf('\n == Benchmark - time LRT_recon : %4.2f == \n\n',toc(benchstart))