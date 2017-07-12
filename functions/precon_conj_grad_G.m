function Gk=precon_conj_grad_G(G,C,A,Y,alpha,Psi,d,Phi,F)
tic; 
fprintf('---conjugate gradient algorithm: solving for G_k ---- \n')

tol=1e-10;
maxiter=50; %temp

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
mfun=@(x) Res(bsxfun(@times, ResA(x),w)); %test

% try to reshape operator so we have proper matrix, vector calculations
X= @(G) (L'*(L*G)) +a2PP*G;
Res= @(x) reshape(x,[numel(G),1]);
ResA= @(x) reshape(x,size(G));
Aop = @(G) Res(X(ResA(G))); 



[x,flag,relres,iter,resvec]=bicgstab(Aop,b,tol,maxiter,mfun,[],G(:)); %add initial guess


Gk=ResA(x);
figure(19); plot(log10(resvec)); 

t=toc; 
fprintf('Time taken: %d seconds',t)
fprintf('| relres %d | iters: %d | \n',relres,iter)

end
