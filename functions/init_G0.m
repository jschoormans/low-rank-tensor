function [G0,C0,A,B,Y,Z]=init_G0(P,Phi,L1)
%initializes all matrices.

%input: 
% P: 1-unfolded zero-filled recon
% Phi: kronecker product of subspace estimates
% L1: estimated rank of G
disp('initializing G_0...')
X=P*Phi'; 
[U,S,V] = svd(X); 
G0=U(:,1:L1);                           % first Lg vectors from left-dominant  svd of P10 Psi^H

C0=G0'*P*Phi';                          % G0^H P10 Psi^H 

A = zeros(size(G));
B = zeros(size(C));
Y = zeros(size(G));
Z = zeros(size(C));
end
