

%% set breakpoint before int G0


%>>>>>>>>>>>>> init G0 part

Phi=kron(nav_estimate_2,nav_estimate_1).';      %from subspaces Phi= kron(G^4,G^3)^T


P0=F'*(kspace);
P1_0=reshape(P0,unfoldedIsize); %1-unfolding of zero filled recon (what is the shape of this matrix?)

X=P1_0*Phi'; 
% [U,S,V] = svds(X,L1);
[U,S,V] = svd(X,'econ'); 

G0=U(:,1:params.Lg);                           % first Lg vectors from left-dominant  svd of P10 Psi^H
% this is the problem isnt it: G0 is actually G0C0 ---> scaling completely
% wrong 

C0=G0'*P1_0*Phi';                          % G0^H P10 Psi^H 
d=reshape(kspace,[430*129*4,6*48]);

%>>>>>>>>>>>>>.........
samplingmask=abs(d)>0; 
L=Lambda(F,C0,Phi,samplingmask);

G_est=L'*d;

figure(1); clf; plot(abs(G_est(:,1))); hold on; plot(abs(G0(:,1)));hold off; 
legend('b1 (L^T*data)','b2 :G0');

% this is an issue!!!

plot(abs(G_est)./abs(G0))
title('G_est / G')

%%%==> something wrong in lambda operator. Every spatial rank is scaled
%%% different. 
% In CG: this means that L'd is overrepresented wrt the contribution of A_k
% 


%% test if problem is in F operator
P1_0=reshape(P0,unfoldedIsize); %1-unfolding of zero filled recon (what is the shape of this matrix?)
P1_0_2=F'*d;
sum(sum(P1_0~=P1_0_2)); %=0: exactly the samae

%% test if Lambda is not ok...
% should be true: L *L'd  ~ d 
% should be true : L'*L*G = G

G2=L'*(L*G);
d2=L*(L'*d);


figure(3); plot(abs(G2(:)./abs(G(:))));
figure(4); plot(abs(d2(:)./abs(d(:))));
%==> lambda not okay at all... 

%%  redo the d example

G1= F'*(d*Phi'*pinv(C))%  GC 
d1=(F*(G1*C*Phi)); %sure


figure(4);clf;  plot(abs(d1(:,288)./abs(d(:,288))));
