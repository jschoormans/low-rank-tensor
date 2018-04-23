function params= params_init(kspace)
% initializes relevant parameters for the LRT reconstruction 
% possible to change these before running recon

params.sizes.kspace = size(kspace);
vel_enc_directions = 1;

params.inspectLg=false;              % option to dynamically choose spatial rank based on first guess of C
params.increase_penalty_parameters=true; %option to increase alpha and beta inside the loop
params.scaleksp=true;              %option to scale kspace before recon - for consistent parameter use (TODO!)

params.Lg=6;                        %spatial rank
params.L3=6;                        %rank of first parameter dimension (cardiac)
params.L4=vel_enc_directions;                        %rank of second parameter dimension (vel.enc.)

params.sparsity_transform='TV';     %'TV'/'wavelet'/'TVOP'/'I'
params.niter=20;                    %number of outer iterations in algo

%initialize parameters
params.alpha= 1;                    %penalty parameter >0
params.beta=  1;                    %penalty parameter >0
params.lambda=0.8483;                 %sparsity parameter
params.mu=0.8483;                      %sparsity parameter

params.Imref=[];                    %possible reference (gold standard) image
params.x=35;                        %pixel to plot during recon loop
params.y=52;                        %pixel to plot during recon loop

params.C.tol=1e-10;
params.C.maxiter=50;
params.G.tol=1e-13;
% params.G.maxiter=8;
params.G.maxiter=50;
% params.G.precon=false;               %optional preconditioning of G
params.G.precon=true;

% params.subspacedim1=1               % dimension along which to take f.s. vals 
% params.subspacedim2=1               % dimension along which to take f.s. vals
params.columns=1:vel_enc_directions;  % max number of columns (extension of subspacedim1 and 2)
params.rows=1:24;                     % max number of rows
% params.columns=2;  % max number of columns (extension of subspacedim1 and 2)
% params.rows=8;                     % max number of rows


params.nav_estimate_1=[];
params.nav_estimate_2=[];
params.eigenvals_1=[];
params.eigenvals_2=[];

params.autolambda=0   ; 
params.automu=0       ;             % automatically estimate mu on s.t. of first iter
params.normalize_sense=1;           %automatically normalizes sense maps 

params.hadamard = 0;                % Whether or not to exploit Hadamard sparsity
params.visualize = 0; % does not work yet in parameter_optimization

params.nullbackground = 0;
