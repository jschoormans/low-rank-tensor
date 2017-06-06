% making undersampled k-space measurement
% with fully sampled navigator k-space centers 
% var density sampling for all other lines 
% function or not?? 
% TO DO: 

rng(3); % random seed to control phantom generation
res=128;

ADCvals=[0.1 0.009 1 3 6];
T1vals=[20 50 100 200 800].*1e-3;

TE=[10:10:100].*0.001 % in seconds
bvals=3*[0:100:900].*1e-3 %units?

I=diffusion_T1_phantom(res,ADCvals,T1vals,TE,bvals,1);

figure(1); Q=[]
for ii=1:size(I,3)
    J=[];
    for jj=1:size(I,4);
        J=[J,abs(I(:,:,ii,jj))];
    end
    Q=[Q;J];
end
imshow(abs(Q),[0 1])
xlabel('echo times')
ylabel('b-values')

%% to k space 
d=fftshift(fftshift(fft(fft(ifftshift(ifftshift(I,2),1),[],1),[],2),1),2);
figure(2); imshow(abs(d(:,:,1,1)),[])

%% make undersampling mask  
uf=0.5; % undersampling factor 
mask=rand(size(d))>(1-uf); %undersampling

% add center
ctr=10;
ctrcoords=floor(res/2)-ctr: floor(res/2)+ctr;
ll=length(ctrcoords);
mask(ctrcoords,ctrcoords,:,:)=ones(ll,ll,size(d,3),size(d,4));

figure(3); imshow(mask(:,:,1,1))

%% make undersampled measurement; 
d_u=mask.*d; 












