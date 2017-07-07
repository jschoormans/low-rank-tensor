
tic; 
fprintf('---conjugate gradient algorithm: solving for G_k ---- \n')

tol=1e-10;

L=Lambda(F,C,Phi,abs(d)>0);
a2PP=(alpha/2)*Psi'*Psi;

b=(L'*d + (alpha/2)*Psi'*(A+Y./alpha));
b=b(:);
% try to reshape operator so we have proper matrix, vector calculations

X= @(G) (L'*(L*G)) +a2PP*G;
Res= @(x) reshape(x,[numel(G),1]);
ResA= @(x) reshape(x,size(G));
Aop = @(G) Res(X(ResA(G))); 

maxiter=200;

% m1fun=@(G) Res(fft2c(reshape(G,[128 128 6]))); %test
% m2fun=@(G) Res(ifft2c(reshape(G,[128 128 6]))); %test
CPPC=(C*Phi*Phi'*C');
iCPPC=inv(CPPC);
w=diag(CPPC)+alpha/2; 
w=permute(w,[2 1]);
% w=[1,2,3,4,5,6].^2
w=1./w;

m1fun=@(x) Res(bsxfun(@times, ResA(x),w)); %test
m2fun=@(x)  Res(bsxfun(@times, ResA(x),1)); %test

[xpcg,flag,relres,iter,resvecpcg]=bicgstab(Aop,b,tol,maxiter,m1fun,[],G(:)); %add initial guess
[x,flag,relres,iter,resvec]=bicgstab(Aop,b,tol,maxiter,[],[],G(:)); %add initial guess

clf
figure(19); plot(log(resvec),'ko'); 
hold on; plot(log(resvecpcg),'r*');
xlabel('iter'); ylabel('log(error)')
legend('normal','pcg')