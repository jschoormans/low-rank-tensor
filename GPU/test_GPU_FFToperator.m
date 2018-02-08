
%%
ll=44
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

temp=[1:ll]; temp=2*(mod(temp,2)-0.5).* 1i;

fprintf('bsxfun \n')
tic
res=bsxfun(@mtimes,bsxfun(@mtimes,x,temp),temp.');
res = fft2(res);
res_bsxfun=bsxfun(@mtimes,bsxfun(@mtimes,res,temp),temp.');
toc


figure(1); clf;hold on; 
plot(abs(res_bsxfun(:,round(ll/2),1,1,1,1)),'k-')
plot(abs(res_fft2(:,round(ll/2),1,1,1,1)),'r+')
legend('bsxfun','fftshift')

%% time bsxfun versus fftshift for a number of sizes



%% check for ifft

for ll=[2:2:80]
A=rand(ll,ll,100,20,5);
x=gpuArray(A); 

fprintf('improved with fft2\n ')
%new improved with fft2 %0.43 s
tic
res_fft2 = fftshift(fftshift(fft2(ifftshift(ifftshift(x,1),2)),1),2);
t_fft2(ll)=toc


temp=[1:ll]; temp=2*(mod(temp,2)-0.5).* 1i; %can we put this outside the iterations?? yes
fprintf('bsxfun \n')
tic
res=bsxfun(@mtimes,bsxfun(@mtimes,x,temp),temp.');
res = fft2(res);
res_bsxfun=bsxfun(@mtimes,bsxfun(@mtimes,res,temp),temp.');
t_bsxfun(ll)=toc
end


%%
figure(2); clf; hold on; plot([2:2:80],t_bsxfun(t_bsxfun>0),'k-+'); plot([2:2:80],t_fft2(t_fft2>0),'r-+');
legend('bsxfun','fftshift')
xlabel('N'); ylabel('time for whole FFT operation')

%% IFFTMOD


ll=40
A=rand(ll,ll,100,20,5);
x=gpuArray(A); 


fprintf('improved with fft2\n ')
%new improved with fft2 %0.43 s
tic
res_fft2 = fftshift(fftshift(ifft2(ifftshift(ifftshift(x,1),2)),1),2);
toc

temp=[1:ll]; temp=2*(mod(temp,2)-0.5).* 1i;
temp=conj(temp); 
fprintf('bsxfun \n')
tic
res=bsxfun(@mtimes,bsxfun(@mtimes,x,temp),temp.');
res = ifft2(res);
res_bsxfun=bsxfun(@mtimes,bsxfun(@mtimes,res,temp),temp.');
toc


figure(3); clf;hold on; 
plot(imag(res_bsxfun(:,round(ll/2),1,1,1,1)),'k-')
plot(imag(res_fft2(:,round(ll/2),1,1,1,1)),'r+')
legend('bsxfun','ifftshift')

sum(res_bsxfun(:)==res_fft2(:))
numel(res_bsxfun)