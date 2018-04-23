function Ck=precon_conj_grad_C(G,C,B,Z,beta,d,Phi,F,params)
tic; 
fprintf('CG for C_k: ')
samplingmask=abs(d)>0; 
tol=params.C.tol;
maxiter=params.C.maxiter;

if params.hadamard==0;
    L=Lambda2(F,G,Phi,samplingmask);
else
    L=Lambda2_hadam(F,G,Phi,samplingmask,params);
end

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
% line 26 to 29 removed for now 15-11-2017, 30-11-2017 added
if params.visualize == 1;
    figure(998);subplot(224);
    plot(log10(resvecpcg./norm(b(:))),'r*-'); 
    xlabel('iterations'); ylabel('10log of residual')
    drawnow; 
end
t=toc; 
fprintf('t: %i seconds',t)
fprintf('| relres %d | iters: %i | \n',relres,iter)

end
