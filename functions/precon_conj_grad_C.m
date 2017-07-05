function Ck=precon_conj_grad_C(G,C,B,Z,beta,d,Phi,F)

tol=1e-19;
maxiter=100;

L=Lambda2(F,G,Phi,abs(d)>0);

b=((L'*d) + (beta/2)*(B+Z./beta));
b=b(:);

X= @(G) (L'*(L*G)) +(beta/2)*G;
Res= @(x) reshape(x,[numel(C),1]);
ResA= @(x) reshape(x,size(C));
Aop = @(C) Res(X(ResA(C))); 

%%
% [x2,flag,relres,iter,resvec2]=bicgstab(Aop,b,tol,maxiter);
% [x2,flag,relres,iter,resvec2]=pcg(Aop,b,tol,maxiter);
[x2,flag,relres,iter,resvec2]=cgs(Aop,b,tol,maxiter); %seems fastest

Ck=ResA(x2);
figure(29); plot(log(resvec2)); drawnow; 


end
