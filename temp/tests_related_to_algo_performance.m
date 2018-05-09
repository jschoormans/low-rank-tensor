
%% check relative signal contributions in CG 
%% set breakpoint in CG algo (after calcualtion of b) 
b1=L'*d;
b2=(alpha/2)*(Psi'*(A+Y./alpha));
b1=gather(b1);
b2=gather(b2); 

figure(1); clf; plot(abs(b1(:,1))); hold on; plot(abs(b2(:,1)));hold off; 
legend('b1 (L^T*data)','b2 :(alpha/2)*Psi^T*(A+Y/alpha)');

%% check outer iteration
% set breakpoint after 2 iters or recon
size(Ak)
size(Psi*G)
size(Y)

PsiG=Psi*G; 
figure(1); clf; 
plot(abs(Ak(:,1)),'-'); hold on; 
plot(abs(PsiG(:,1)),'.r');
plot(abs(Y(:,1)./ alpha),'+b');hold off; 
legend('Ak','PsiG','Y/alpha'); 

%
L=Lambda(F,C,Phi,mask);
d=reshape(kspace,[430*129*4,6*48]);
G_est=L'*d;

figure(5); clf;
plot(abs(G(:,1)),'.r'); hold on
plot(abs(G_est(:,1)),'.b'); hold off
legend('G in algo','G from kspace')

figure(6); clf;
plot(abs(G_est(:))./abs(G(:)),'-r');
title('G_est divided by G')


%% check Lambda operator 
%% set breakpoint after estimation of G 

L=Lambda(F,C,Phi,mask);
d=reshape(kspace,[430*129*4,6*48]);
G_est=L'*d;

size(G_est)
size(G)

figure(1); clf; plot(abs(G_est(:,1))); hold on; plot(abs(G(:,1)));hold off; 
legend('b1 (L^T*data)','b2 :G');

% this is an issue!!!

plot(abs(G_est)./abs(G))
title('G_est / G')

%%%==> something wrong in lambda operator. Every spatial rank is scaled
%%% different. 
% In CG: this means that L'd is overrepresented wrt the contribution of A_k
% 


