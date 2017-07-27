function Gk=precon_conj_grad_G(G,C,A,Y,alpha,Psi,d,Phi,F,params)
tic; 
fprintf('CG for G_k: ')

tol=params.G.tol;
maxiter=params.G.maxiter;

L=Lambda(F,C,Phi,abs(d)>0);
a2PP=(alpha/2)*Psi'*Psi;

%input data
b=(L'*d + (alpha/2)*Psi'*(A+Y./alpha));
b=b(:);

%build preconditioner
CPPC=(C*Phi*Phi'*C');
w=diag(CPPC)+alpha/2; 
w=permute(w,[2 1]);
w=1./w;

% try to reshape operator so we have proper matrix, vector calculations
X= @(G) (L'*(L*G)) +a2PP*G;
Res= @(x) reshape(x,[numel(G),1]);
ResA= @(x) reshape(x,size(G));
Aop = @(G) Res(X(ResA(G))); 
mfun=@(x) Res(bsxfun(@times, ResA(x),w)); %test

[x,flag,relres,iter,resvec]=bicgstab(Aop,b,tol,maxiter,mfun,[],G(:)); %add initial guess

Gk=ResA(x);
figure(19); plot(log10(resvec./norm(b(:))),'r*-'); 
xlabel('iterations'); ylabel('10log of relative residual')
t=toc; 
fprintf('Time taken: %d seconds',t)
fprintf('| relres %d | iters: %d | \n',relres,iter)

end
