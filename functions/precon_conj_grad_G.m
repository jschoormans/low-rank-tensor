function Gk=precon_conj_grad_G(G,C,A,Y,alpha,Psi,d,Phi,F)
% TO DO: THINK VERY WELL ABOUT MASK (OMEGA) 
tol=10;
maxiter=10;

L=Lambda(F,C,Phi,abs(d)>0);
a2PP=(alpha/2)*Psi'*Psi;

b=(L'*d + (alpha/2)*Psi'*(A+Y./alpha));
b=b(:);
% try to reshape operator so we have proper matrix, vector calculations

X= @(G) (L'*(L*G)) +a2PP*G;
Res= @(G) reshape(G,[102400*24,1]);
ResA= @(G) reshape(G,[102400,24]);
Aop = @(G) Res(X(ResA(G))); 


[x2,flag,relres,iter,resvec2]=bicgstab(Aop,b);
% [x1,flag,relres,iter,resvec1]=pcg(Aop,b,[],20);

Gk=ResA(x2);
figure(19); plot(resvec2); drawnow; 

%% TO DO:calculate preconditioner
% CAN BE ADDED AS FUNCTION 

%{
x=G(:);
R=b(:)-Aop(x);
P=R(:);


iter=0;resid=1e99;
while iter<maxiter && resid>tol
    iter=iter+1;
    alpha= (R'*R)/(P'*Aop(P));
    xk=x+alpha.*P;
    Rk=R-alpha.*(Aop(P));
    
    iters=    [0:0.2:2]
    for ii=1:length(iters)
        
    nn(ii)=norm(Aop(x+(alpha.*iters(ii)).*P)-b(:),2);
    end
    figure(19); plot(iters,nn); drawnow
    
    resid=abs(sum(Rk(:)));
    disp(['iteration: ',num2str(iter),'|residual: ',num2str(resid)])
    
    beta=(Rk'*Rk)/(R'*R);
    Pk=Rk+beta.*P;
    R=Rk; P=Pk; x=xk;
end
Gk=xk;
%}
end
