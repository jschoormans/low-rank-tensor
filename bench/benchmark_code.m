% recon scan 16 (high res)
tic
fprintf('\n == Benchmark == \n',toc)
fprintf('Loading data... \n',toc)

load('benchmark-data.mat')
% load
params.visualization=1; 
params= params_init(kspaceinput)
params.GPU=1;
params.rows=[1:6]
params.columns=[1:48]
params.TVoption=4 %for CPU

params.C.maxiter=50
params.G.maxiter=50

fprintf('Starting recon... \n',toc)
benchstart=tic;
params.GPUdouble=1;         %for GPUs supporting double (FLUX!!)
P_recon=LRT_recon(kspaceinput,squeeze(sens),params);

fprintf('\n == Benchmark - time LRT_recon : %4.2f == \n\n',toc(benchstart))