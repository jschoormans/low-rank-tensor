% recon scan 16 (high res)
tic
fprintf('\n == Benchmark == \n',toc)
fprintf('Loading data... \n',toc)

load('benchmark-data.mat')
% load
params.GPU=1;
params.visualization=1; 
params.niter=15
params.x=28; params.y=125;
params.increase_penalty_parameters=1
params.Lg=12;
params.lambda=0.5
params.alpha=1e2
params.mu=2

fprintf('Starting recon... \n',toc)
benchstart=tic;
params.GPUdouble=0; %for GPUs supporting double 
[P_recon,params_out]=LRT_recon(kspaceinput,squeeze(sens),params);

fprintf('\n == Benchmark - time LRT_recon : %4.2f == \n\n',toc(benchstart))

imagine(abs(squeeze(P_recon(:,:,:,:,5))))
%% project recon onto subspace (project on dim2) - SET BREAKPOINT AT P0 in LRT_Recon

if 0
    Itemp=squeeze(P0(:,:,1,2,:));
    A=params.nav_estimate_2; %subspace
else
    Itemp=squeeze(P0(:,:,1,:,24));
    A=params.nav_estimate_1; %subspace
end
Itemp=reshape(Itemp,[size(Itemp,1)*size(Itemp,2),size(Itemp,3)]); %make into matrix 

size(Itemp)
size(A)

Pr=Itemp*A;
size(Pr)

Prr=reshape(Pr,[size(P0,1),size(P0,2)*size(Pr,2)]);
Prr=Prr./max(Prr(:));

figure(1);
imshow(abs(Prr(:,:)),[])

