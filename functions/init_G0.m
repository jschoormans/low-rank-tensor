function G0=init_G0(P,Phi,L1)
%input: 
% P: 1-unfolded zero-filled recon
% Phi: kronecker product of subspace estimates
% L1: estimated rank of G
disp('initializing G_0...')
X=P*(Phi); 
[U,S,V] = svd(X); %econ because 
G0=U(:,1:L1);

end