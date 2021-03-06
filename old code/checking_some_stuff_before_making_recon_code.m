% try  diffusion and T1 phantom 
% makes image for different b values and echo times
% then this image tensor is described as a low-rank tensor 

% phantom parameters
res=128;

ADCvals=[0.1 0.009 1 3 6];
T1vals=[20 50 100 200 800].*1e-3;

TE=[10:10:100].*0.001 % in seconds
bvals=3*[0:100:900].*1e-3 %units?

I=diffusion_T1_phantom(res,ADCvals,T1vals,TE,bvals,1);
%
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

%% kspace : called d
d=fftshift(fftshift(fft(fft(ifftshift(ifftshift(I,2),1),[],1),[],2),1),2);
figure(2); imshow(abs(d(:,:,1,1)),[])
figure(11); Q=[]
for ii=1:size(I,3)
    J=[];
    for jj=1:size(I,4);
        J=[J,abs(d(:,:,ii,jj))];
    end
    Q=[Q;J];
end
imshow(abs(Q),[0 1])
xlabel('echo times')
ylabel('b-values')
%% TUCKER DECOMPOSITION

rank1=60
rank2=60
rank3=6
rank4=6
[F1,C,pctage,errormap]=tucker(I,[rank1,rank2,rank3,rank4]);


% check d-omega(F G C Psi)
% omega is selection operator???
% psi = kronecker  product of G^{n}, G^(n-1)... (all but the spatial dimensions)

Psi=kron(F1{4},F1{3}) %kronecker product of G^3 and G^4 (echo times and b-values)
G=kron(F1{2},F1{1});  % spatial low-rank matrix 
C1=reshape(C,[rank1*rank2,rank3*rank4]); %flattened C 

image_recon=(G*C1*Psi.');
image_tensor_res=reshape(image_recon,[res,res,10,10]);
kspace_recon=fftshift(fftshift(fft(fft(ifftshift(ifftshift(image_tensor_res,2),1),[],1),[],2),1),2);

%test image_recon
figure(2); 
subplot(121); imshow(image_tensor_res(:,:,1,1),[])
subplot(122); imshow(I(:,:,1,1),[])

l2_difference=sqrt(sum(abs(d(:)-kspace_recon(:)).^2))
l2_only_d=sqrt(sum(abs(d(:)).^2))
l2_difference/l2_only_d

%% CAN WE MAKE PSI FROM THE NAVIGATORS??

nav_parameter_dim1 = squeeze(d(54:74,54:74,:,5));
nav_parameter_dim2 = squeeze(d(54:74,54:74,5,:));

figure(10); subplot(211);
montage(permute(abs(nav_parameter_dim1),[1 2  4 3]),'displayrange',[],'size',[2 5]);
subplot(212);
montage(permute(abs(nav_parameter_dim2),[1 2  4 3]),'displayrange',[],'size',[2 5]);

S1=reshape(nav_parameter_dim1,[21*21,10]);
S2=reshape(nav_parameter_dim2,[21*21,10]);

[left_1,eigen_1,right_1]=svd(S1);

nav_estimate_1=right_1(:,1:rank3);

figure(11); subplot(121); imshow(F1{3},[]);
subplot(122); imshow(nav_estimate_1,[]);
theta = subspace(F1{3},nav_estimate_1) %angle between two subspaces 
figure(12); hold on; 
plot(abs(F1{3}(:,1)),'-*g'); 
plot(abs(nav_estimate_1(:,1)),'-*b');
hold off

%% make PSI from navigator function
nav_estimate= subspace_estimator(nav_parameter_dim1,4);
figure(99); imshow(nav_estimate);

