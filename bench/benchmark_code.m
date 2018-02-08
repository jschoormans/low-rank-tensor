% recon scan 16 (high res) 
fprintf('\n == Benchmark == \n',toc)

load('benchmark-data.mat')
% load
params.GPU=0;

benchstart=tic;
P_recon=LRT_recon(kspaceinput(:,1:128,:,:,:),squeeze(sens(:,:,1:128,:)),params);

fprintf('\n == Benchmark - time LRT_recon : %4.2f == \n\n',toc(benchstart))