
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
m3fun=@(x) D2inv*x %(run other code first)

[xpcg,flag,relres,iter,resvecpcg]=bicgstab(Aop,b,tol,maxiter,m1fun,[],G(:)); %add initial guess
[x,flag,relres,iter,resvec]=bicgstab(Aop,b,tol,maxiter,[],[],G(:)); %add initial guess
[xpcg2,flag,relres,iter,resvecpcg2]=bicgstab(Aop,b,tol,maxiter,m3fun,[],G(:)); %add initial guess

clf
figure(19); plot(log10(resvec),'ko'); 
hold on; plot(log10(resvecpcg),'r*');
plot(log10(resvecpcg2),'g*'); 
xlabel('iter'); ylabel('log(error)')
legend('normal','Jacobi preconditioning','full inverse preconditioning');
hold off


%% TEST advanced preconditioning 12 - 7 -2017
n2=size(G,1) %number of pixels 
Gvec=Res(G);
size(Gvec) % size: spatial rank * npixels 
size(CPPC)

%so F*FGCPPC + a/lpa/2 ~= GCPPC + alpha/2

% define pseudeo operator D

D11=opDiag(repmat(CPPC(1,1),[1,n2]))
D12=opDiag(repmat(CPPC(1,2),[1,n2]))
D13=opDiag(repmat(CPPC(1,3),[1,n2]))
D21=opDiag(repmat(CPPC(2,1),[1,n2]))
D22=opDiag(repmat(CPPC(2,2),[1,n2]))
D23=opDiag(repmat(CPPC(2,3),[1,n2]))
D31=opDiag(repmat(CPPC(3,1),[1,n2]))
D32=opDiag(repmat(CPPC(3,2),[1,n2]))
D33=opDiag(repmat(CPPC(3,3),[1,n2]))
D=[D11,D12,D13;D21,D22,D23;D31,D32,D33]
alphaop=opDiag(repmat(alpha/2,size(Gvec)));
D2=D+alphaop

resultD=D*Gvec;
resultD2=D2*Gvec;
result_compare=Aop(G);
figure(19); clf; hold on; plot(abs(result_compare)); plot(abs(resultD)); plot(abs(resultD2)); hold off
legend('A operator','D without alpha','D with alpha')


D2inv=opInverse(D2);
resultD2inv=D2inv*Gvec;
resultD2inv=ResA(resultD2inv);
figure(21)
imshow(reshape(resultD2inv(:,1),[128 128]),[]);
title('inverse of pseudo-operator matrix')

%% 
%LU decomposition










