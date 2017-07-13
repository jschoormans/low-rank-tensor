% making undersampled k-space measurement
% with fully sampled navigator k-space centers 

% TO DO: -add var density and more us options
rng('default') 
rng(3); % random seed to control phantom generation
res=128;


%>>>>>>>>>>>>>>> generate Diffusion & T1 (T2?) weigthed phantoms<<<<<<<<<<<
% ADCvals=[1:5].*1e-3;
% T1vals=[20 50 100 200 800].*1e-3;
% 
% TE=[10:10:100].*0.001; % in seconds
% bvals=[0:100:900]; %units?
% 
% I=diffusion_T1_phantom(res,ADCvals,T1vals,TE,bvals,1);


%>>>>>>>>>>>>>> generate VFA TSE & T2 weighted phantoms<<<<<<<<<<<<<<<<<<<<
T1vals=[500 400 300 200 1000 500 400 300 200 1000].*1e-3; % in seconds
T2vals=[20 30 40 50 50 55 50 60 70 10].*1e-3; % in seconds

T2prep=[10:10:97].*0.001; % in seconds
TEes=4.*0.001; % in seconds
ETL = 20; 

I=VFA_TSE_T2_T1_phantom(res,T2vals,T1vals,T2prep,TEes,ETL,1,complexsim);


figure(1); Q=[];
for ii=1:size(I,3)
    J=[];
    for jj=1:size(I,4);
        J=[J,abs(I(:,:,ii,jj))];
    end
    Q=[Q;J];
end
imshow(abs(Q),[0 1])
xlabel('Prameter 1')
ylabel('Prameter 2')
clear Q J 
%%  add coil sensitivity information 
[Ic,sens]=addcoilsensitvity_to_simulated(I,ncoils);

%% to k space 
d=fftshift(fftshift(fft(fft(ifftshift(ifftshift(Ic,2),1),[],1),[],2),1),2)./sqrt(res*res);
figure(2); imshow(abs(d(:,:,1,1)),[]); axis off; title('k-space')

%% make undersampling mask  (independent of coil dimensions!)
mask=rand(size(I))>(1-uf); %undersampling

% add center for subspace estimae
ctr=10;
ctrcoords=floor(res/2)-ctr: floor(res/2)+ctr; 
ll=length(ctrcoords);
mask(ctrcoords,ctrcoords,1,:)=ones(ll,ll,1,size(I,4));
mask(ctrcoords,ctrcoords,:,1)=ones(ll,ll,size(I,3),1);

% for all measurements
ctrcoordsall=floor(res/2)-3: floor(res/2)+3; 
mask(ctrcoordsall,ctrcoordsall,:,:)=ones(7,7,size(I,3),size(I,4));

figure(3); Q=[];
for ii=1:size(I,3)
    J=[];
    for jj=1:size(I,4)
        J=[J,abs(mask(:,:,ii,jj))];
    end
    Q=[Q;J];
end
imshow(abs(Q),[0 1]);axis off;  clear Q J 

mask=repmat(mask,[1 1 1 1 ncoils]);
mask=permute(mask,[1 2 5 3 4]);
%% make undersampled measurement; 
du=mask.*d; 
clear d 

%add noise
if noiselevel>0
du=du+(randn(size(du)).*mean(du(:)).*noiselevel).*(du~=0);
end









