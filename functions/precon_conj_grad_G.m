function Gk=precon_conj_grad_G(G,C,A,Y,alpha,Psi,d,Phi,F)
% TO DO: THINK VERY WELL ABOUT MASK (OMEGA) 
tol=1e-13;
maxiter=200;

L=Lambda(F,C,Phi,abs(d)>0);
a2PP=(alpha/2)*Psi'*Psi;

b=(L'*d + (alpha/2)*Psi'*(A+Y./alpha));
b=b(:);
% try to reshape operator so we have proper matrix, vector calculations

X= @(G) (L'*(L*G)) +a2PP*G;
Res= @(x) reshape(x,[numel(G),1]);
ResA= @(x) reshape(x,size(G));
Aop = @(G) Res(X(ResA(G))); 

[x2,flag,relres,iter,resvec1]=cgs(Aop,b,tol,maxiter); %seems fastest
% [x2,flag,relres,iter,resvec2]=bicgstab(Aop,b,tol,maxiter);
% [x2,flag,relres,iter,resvec3]=pcg(Aop,b,tol,maxiter);

Gk=ResA(x2);
figure(19); plot(log(resvec1)); 
drawnow;
%% TO DO:calculate preconditioner
% CAN BE ADDED AS FUNCTION 

end
