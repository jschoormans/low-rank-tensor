% recon scan 16 (high res)
vars;
addpath(genpath('L:\basic\divi\Projects\cosart\tensor'));
addpath(genpath('C:\Users\jschoormans\gpuNUFFT-master'))
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\imagine'));


cd('L:\basic\divi\Projects\cosart\tensor\low-rank-tensor\radial')
load synthetic_radial.mat



% load

a=gpuArray(2); clear a; %stupid matalab error workaround 

osf =1; wg = 3; sw = 8;
MDFT = MDgpuNUFFT(trajGPU,w,osf,wg,sw,[res,res,1],[],true);

radialksp= MDFT*I_sorted;
tempI=MDFT'*radialksp;
%%

params=params_init();
params.GPU=1;
params.visualization=1; 
params.TVoption=4;
params.GPUdouble=0; %for GPUs supporting double 
params.nx=128;

% SVD stuff (works better with multicoil data
for ii=1:size(radialksp,3)
    for jj=1:size(radialksp,4)
        coords_center=[64:128:128*11]
        centerksp(:,1,ii,jj)=radialksp(coords_center,1,ii,jj);
    end
end

matrix_dim3=reshape(permute(centerksp,[1 2 4 3]),[11*5,10]);
matrix_dim4=reshape(centerksp,[11*10,5]);
[nav_estimate_1,params.eigenvals_1]= subspace_estimator_multicoil(matrix_dim3,params.L3);

[nav_estimate_2,params.eigenvals_2]= subspace_estimator_multicoil(matrix_dim4,params.L4);

params.nav_estimate_1=nav_estimate_1;
params.nav_estimate_2=nav_estimate_2;
params.params.eigenvals_1=params.eigenvals_1;
params.params.eigenvals_2=params.eigenvals_2;




%%
%TO DO: ADD weight in operators (sqrt(w)) + test parfor...


params.GPU=1;
params.visualization=1; 
params.TVoption=2;

params.visualization=1
params.Lg=6
params.L4=3
params.L3=3
params.G.maxiter=3
params.inspectLg=0
params.niter=12;
params.G.precon=0
params.G.maxiter=10
params.C.maxiter=10
params.lambda=0.05;
params.alpha=1.2
params.mu=1.2

P_recon=LRT_recon_radial(radialksp,[],params,MDFT);
