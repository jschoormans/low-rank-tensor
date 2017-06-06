function G0=init_G0(P,Psi,L1)
%input: 
% P: 1-unfolded zero-filled recon
% Psi: kronecker product of subspace estimates
% L1: estimated rank of G
disp('initializing G_0...')
X=P*conj(Psi); 
[U,S,V] = svd(X); %econ because 
G0=U(:,1:L1);

end
