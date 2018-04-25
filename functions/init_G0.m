function [Phi,G0,C0,A,B,Y,Z]=init_G0(P,Psi,params)
%initializes all matrices.
nav_estimate_1 = params.nav_estimate_1;
nav_estimate_2 = params.nav_estimate_2;
L1 = params.Lg;
%input: 
% P: 1-unfolded zero-filled recon
% Phi: kronecker product of subspace estimates
% L1: estimated rank of G
disp('initializing matrices...')

if params.indepvenccols;
    Phi=nav_estimate_1.';
else
    Phi=kron(nav_estimate_2,nav_estimate_1).';      %from subspaces Phi= kron(G^4,G^3)^T
end

X=P*Phi'; 
% [U,S,V] = svds(X,L1);
[U,S,V] = svd(X,'econ'); 
%Added by Bobby 24-01-2018
eigenvals_initG=diag(S);  %output for evaluation
figure(1234); hold on; plot(eigenvals_initG./max(eigenvals_initG(:)),'r'); plot(L1,eigenvals_initG(L1)./max(eigenvals_initG(:)),'ro');hold off;
title('eigenvalues for the spatial matrix');

G0=U(:,1:L1);                           % first Lg vectors from left-dominant  svd of P10 Psi^H

C0=G0'*P*Phi';                          % G0^H P10 Psi^H 

A=zeros(size(Psi*G0));
Y=zeros(size(A));

B = zeros(size(C0));
Z = zeros(size(C0));
end
