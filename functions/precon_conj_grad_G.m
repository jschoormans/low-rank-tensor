function Gk=precon_conj_grad_G(G,C,A,Y,alpha,Psi,d,Phi,F)
tic; 
fprintf('---conjugate gradient algorithm: solving for G_k ---- \n')



tol=1e-10;
maxiter=200;
maxiter=50; %temp

L=Lambda(F,C,Phi,abs(d)>0);
a2PP=(alpha/2)*Psi'*Psi;

b=(L'*d + (alpha/2)*Psi'*(A+Y./alpha));
b=b(:);
% try to reshape operator so we have proper matrix, vector calculations

X= @(G) (L'*(L*G)) +a2PP*G;
Res= @(x) reshape(x,[numel(G),1]);
ResA= @(x) reshape(x,size(G));
Aop = @(G) Res(X(ResA(G))); 


 c [xpcg,flag,relres,iter,resvecpcg]=bicgstab(Aop,b,tol,maxiter,m1fun,m2fun,G(:)); %add initial guess

[x,flag,relres,iter,resvec]=bicgstab(Aop,b,tol,maxiter,[],[],G(:)); %add initial guess

% [x2,flag,relres,iter,resvec1]=cgs(Aop,b,tol,maxiter); %seems fastest
% [x2,flag,relres,iter,resvec2]=bicgstab(Aop,b,tol,maxiter);
% [x2,flag,relres,iter,resvec3]=pcg(Aop,b,tol,maxiter);

Gk=ResA(x);
figure(19); plot(log(resvec)); 
hold on; plot(resvecpcg)

t=toc; 
fprintf('Time taken: %d seconds',t)
fprintf('| relres %d | iters: %d | \n',relres,iter)

end
