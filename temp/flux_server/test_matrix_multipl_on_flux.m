% FLUX SERVER- OUTSTANDING ISSUES 


%% SINGLE/DOUBLE GPUARRAYS

% USIGN DOUBLE - WORKS ON TESLA P100
A=gpuArray(rand([55470,6],'double'));
B=gpuArray(rand([6,288],'double'));

C=A*B;

%% 
%%DOES NOT WORK ON TESLA P100 - DOES WORK ON TITAN XP
% ERROR : 
% % Error using  * 
% % Call to sgemm in CUBLAS failed with error status:
% % CUBLAS_STATUS_EXECUTION_FAILED.
 
A=gpuArray(rand([55470,6],'single'));
B=gpuArray(rand([6,288],'single'));

C=A*B;

%% DOES NOT WORK ON P100 (DOUBLE/SINGLE COMBINATIONS) 
A=gpuArray(rand([55470,6],'double'));
B=gpuArray(rand([6,288],'single'));

C=A*B;

A=gpuArray(rand([55470,6],'single'));
B=gpuArray(rand([6,288],'double'));

C=A*B;

%% BART TOOLBOX TEST - CS RECON

ph=bart('phantom -k -s8 -x256');
poiss=bart('poisson -y1.3 -z1.3 -Y256 -Z256 -C12 -v');
ph=permute(ph,[3 1 2 4]); 
ph=bsxfun(@times,ph,poiss); 
s=bart('phantom -S8 -x256');
size(ph)
size(s)
tic
recon=bart('pics -d5',ph,s);
tCPU=toc; 
size(recon)

figure(1); clf;
subplot(131); imshow(abs(squeeze(poiss)),[]); title('mask')
subplot(132); imshow(abs(squeeze(recon)),[]); title(['pics CPU - t:',num2str(tCPU),'s'])

% gpu recon
tic()
reconGPU=bart('pics -g',single(ph),single(s));
tGPU=toc
subplot(133); imshow(abs(squeeze(reconGPU)),[]); title(['pics GPU - t:',num2str(tGPU),'s'])

%% BART TOOLBOX TEST - SENSITVITY MAPS ESTIMATION 
s=bart('ecalib',ph);        % works 
s=bart('ecalib -a',ph):     % does not work 









