% test arrayfun 
A=phantom(1000); B=gpuArray(A)
%%
% GPU implementation 1 of TV in one direction  (GPUSHIFT function) 

iters=200

tic 
for i=1:iters
B_shift=gpuShift(B,[-1 0]);
TVy=B-B_shift; end
tGPU1=toc(); 

figure(1); subplot(222)
imshow(abs(TVy))

% GPU implementation 2 of TV in one direction  (normal circshift)

tic 
for i=1:iters
B_shift=circshift(B,[-1 0]);
TVy=B-B_shift;
end
tGPU2=toc(); 

figure(1); subplot(223)
imshow(abs(TVy))

% GPU implementation 3 of TV in one direction  GPUshift+ arrayfun

tic 
for i=1:iters
B_shift=gpuShift(B,[-1 0]);
TVy=arrayfun(@minus,B,B_shift); end
tGPU3=toc(); 

figure(1); subplot(224)
imshow(abs(TVy))

% GPU implementation 3 of TV in one direction  normal circshift + arrayfun

tic 
for i=1:iters
B_shift=circshift(B,[-1 0]);
TVy=arrayfun(@minus,B,B_shift); end
tGPU4=toc(); 

% CPU implementation 1 of TV in one direction 
tic
for i=1:iters
A_shift=circshift(A,[-1 0]); 
TVy=A-A_shift; end
tCPU1=toc;

figure(1); subplot(221);  
imshow(abs(TVy))

fprintf('Time CPU 1: %4.8f \n',tCPU1)
fprintf('Time GPU 1: %4.8f - %4.2f times faster \n',tGPU1, tCPU1/tGPU1)
fprintf('Time GPU 2: %4.8f - %4.2f times faster \n',tGPU2,tCPU1/tGPU2)
fprintf('Time GPU 3: %4.8f - %4.2f times faster \n',tGPU3,tCPU1/tGPU3)
fprintf('Time GPU 3: %4.8f - %4.2f times faster \n',tGPU4,tCPU1/tGPU4)


%% testing multicoil TV operator 


ncoils=100
A=phantom(1000);
A=repmat(A,[1 1 ncoils]); %10 coil image

clear tCPU1 tGPU1 tGPU2
for ncoils=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

% 1 for loop  CPU
tic
for ii=1:ncoils
    TVy(:,:,ii)=A(:,:,ii)-circshift(A(:,:,ii),[1 0]); 
end
tCPU1(ncoils)=toc

B=gpuArray(A);
 

tic
for ii=1:ncoils
    TVyGPU1(:,:,ii)=B(:,:,ii)-gpuShift(B(:,:,ii),[1 0]); 
end
tGPU1(ncoils)=toc

tic
for ii=1:ncoils
    TVyGPU1(:,:,ii)=B(:,:,ii)-circshift(B(:,:,ii),[1 0]); 
end
tGPU0(ncoils)=toc
 

tic

TVyGPU2=pagefun(@min,B,gpuShift(B,[1 0]));

tGPU2(ncoils)=toc

end

%%
figure(10);clf; hold on;
plot(tCPU1./tGPU0+eps,'b*-')
plot(tCPU1./tGPU1+eps,'k*-'); plot(tCPU1./tGPU2+eps,'r*-')
xlabel('ncoils')
ylabel('speedup compared to CPU loop')
legend('GPU - circshift & loop','GPU gpuSHIFT & loop','GPU gpuSHIFT & pagefun')
