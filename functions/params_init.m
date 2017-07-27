function params= params_init()
% initializes relevant parameters for the LRT reconstruction 
% possible to change these before running recon

params.Lg=1;                %spatial rank
params.L3=3;                %rank of first parameter dimension
params.L4=3;                %rank of second parameter dimension

params.sparsity_transform='TV';     %'TV'/'wavelet'/'TVOP'
params.niter=20;            %number of outer iterations in algo

%initialize parameters
params.alpha= 2;           %penalty parameter >0
params.beta=  2;           %penalty parameter >0
params.lambda=1e-2;        %sparsity parameter
params.mu=1e1;            %sparsity parameter

params.Imref=[];            %possible reference (gold standard) image
params.x=20;                %pixel to plot during recon loop
params.y=20;                %pixel to plot during recon loop

params.scaleksp=true';      %option to scale kspace before recon - for consistent parameter use (TODO!)

params.C.tol=1e-10;
params.C.maxiter=50;
params.G.tol=1e-13;
params.G.maxiter=100;
