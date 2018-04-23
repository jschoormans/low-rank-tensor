
%%
%WORKS ON TESLA P100
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