function params= params_init()
% initializes relevant parameters for the LRT reconstruction 

params.Lg=1
params.L3=3
params.L4=3
params.sparsity_transform='TV'
params.niter=20  %number of outer iterations in algo

%initialize parameters
params.alpha= 2;           %penalty parameter >0
params.beta=  2;           %penalty parameter >0
params.lambda=1e-8;        %sparsity parameter
params.mu=1e-1;            %sparsity parameter
