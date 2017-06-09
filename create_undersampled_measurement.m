% making undersampled k-space measurement
% with fully sampled navigator k-space centers 

% TO DO: -add var density and more us options

rng(3); % random seed to control phantom generation
res=128;

ADCvals=[1:5].*1e-3;
T1vals=[20 50 100 200 800].*1e-3;

TE=[10:10:100].*0.001 % in seconds
bvals=[0:100:900] %units?

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
clear Q J 
%% to k space 
d=fftshift(fftshift(fft(fft(ifftshift(ifftshift(I,2),1),[],1),[],2),1),2)./sqrt(res*res);
figure(2); imshow(abs(d(:,:,1,1)),[])

%% make undersampling mask  
mask=rand(size(d))>(1-uf); %undersampling

% add center for subspace estimae
ctr=10;
ctrcoords=floor(res/2)-ctr: floor(res/2)+ctr; 
ll=length(ctrcoords);
mask(ctrcoords,ctrcoords,1,:)=ones(ll,ll,1,size(d,4));
mask(ctrcoords,ctrcoords,:,1)=ones(ll,ll,1,size(d,4));

% for all measurements
ctrcoordsall=floor(res/2)-3: floor(res/2)+3; 
mask(ctrcoordsall,ctrcoordsall,:,:)=ones(7,7,size(d,3),size(d,4));

figure(3); Q=[]
for ii=1:size(I,3)
    J=[];
    for jj=1:size(I,4);
        J=[J,abs(mask(:,:,ii,jj))];
    end
    Q=[Q;J];
end
imshow(abs(Q),[0 1]); clear Q J 
%% make undersampled measurement; 
du=mask.*d; 

clear d 

%add noise
if noiselevel>0
du=du+(randn(size(du)).*mean(du(:)).*noiselevel).*(du~=0);
end









