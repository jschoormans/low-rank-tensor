ll=100
A=rand(ll,ll,100,20,5);
A=gpuArray(A); 
%%

tic
F1=fft(fft(A,[],1),[],2);
toc

tic
F2=fft2(A); 
toc


tic
F2=fft2(A); 
toc
%%
ll=98
A=rand(ll,ll,100,20,5);
x=gpuArray(A); 

fprintf('\n \n original method - cpu\n')
tic
res = fftshift(fftshift(fft(fft(ifftshift(ifftshift(A,1),2),[],1),[],2),1),2);
toc

%original method %0.70 s
fprintf('original method - gpu\n')
tic
res = fftshift(fftshift(fft(fft(ifftshift(ifftshift(x,1),2),[],1),[],2),1),2);
toc

fprintf('improved with fft2\n ')
%new improved with fft2 %0.43 s
tic
res_fft2 = fftshift(fftshift(fft2(ifftshift(ifftshift(x,1),2)),1),2);
toc

fprintf('only fft2 - no fftshifts\n')
%only fft2 is %0.14 seconds
tic
res = fft2(x);
toc

temp=[1:ll]; temp=mod(temp,2).* 1i;

fprintf('bsxfun \n')
tic
res=bsxfun(@mtimes,bsxfun(@mtimes,x,temp),temp.');
res = fft2(res);
res_bsxfun=bsxfun(@mtimes,bsxfun(@mtimes,res,temp),temp.');
toc


%%
figure; hold on; 
plot(abs(res_bsxfun(:,1,1,1,1,1)),'k')
plot(abs(res_fft2(:,1,1,1,1,1)),'r')


