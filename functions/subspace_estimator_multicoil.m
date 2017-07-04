function nav_estimate= subspace_estimator_multicoil(kspace,L)
% Function that calculates the subspace based on k-space
% estimator 
% to do: make compatible with 3D kspace 


if ndims(kspace)==4 % kx ky nc param

elseif ndims(kspace)==3 %kx ky param
    kspace=permute(kspace,[1 2 4 3]);    
end  
    
disp('estimating subspace...')
%input:
%kspace: 3D kspace center size(kx,ky,nc,M)
%L:  rank of subspace estimator (L) (L<=M)
assert(L<=size(kspace,4))

%reshape to Casorati matrix
S=reshape(kspace,[size(kspace,1)*size(kspace,2)*size(kspace,3),size(kspace,4)]);

% calculate singular value decomposition
[left_1,eigen_1,right_1]=svd(S);

%navigator estimate are first L left singular vectors
nav_estimate=right_1(:,1:L);



end
