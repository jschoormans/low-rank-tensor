function Ck=precon_conj_grad_C(G,C,B,Z,beta,d,Phi,F)
tic; 
fprintf('---conjugate gradient algorithm: solving for C_k ---- \n')

tol=1e-19;
maxiter=50;

L=Lambda2(F,G,Phi,abs(d)>0);

b=((L'*d) + (beta/2)*(B+Z./beta));
b=b(:);

X= @(C) (L'*(L*C)) +(beta/2)*C;
Res= @(x) reshape(x,[numel(C),1]);
ResA= @(x) reshape(x,size(C));
Aop = @(C) Res(X(ResA(C))); 

%%
% [x2,flag,relres,iter,resvec2]=bicgstab(Aop,b,tol,maxiter);
% [x2,flag,relres,iter,resvec2]=pcg(Aop,b,tol,maxiter);
% [x2,flag,relres,iter,resvec2]=cgs(Aop,b,tol,maxiter); %seems fastest
[xpcg,flag,relres,iter,resvecpcg]=cgs(Aop,b,tol,maxiter,[],[],C(:)); %seems fastest

Ck=ResA(xpcg);
figure(29); plot(log10(resvecpcg./norm(b(:))),'r*-'); 
xlabel('iterations'); ylabel('10log of residual')
drawnow; 

t=toc; 
fprintf('Time taken: %d seconds',t)
fprintf('| relres %d | iters: %d | \n',relres,iter)

end
