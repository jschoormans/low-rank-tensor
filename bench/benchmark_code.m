% recon scan 16 (high res)
tic
fprintf('\n == Benchmark == \n',toc)
fprintf('Loading data... \n',toc)

load('benchmark-data.mat')
% load
params.GPU=1;
params.visualization=1; 
fprintf('Starting recon... \n',toc)
benchstart=tic;
params.GPUdouble=0; %for GPUs supporting double 
P_recon=LRT_recon(kspaceinput,squeeze(sens),params);

fprintf('\n == Benchmark - time LRT_recon : %4.2f == \n\n',toc(benchstart))