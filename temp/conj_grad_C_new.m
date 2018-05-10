function Ck=conj_grad_C_new(G,C,B,Z,beta,d,Phi,F,params)
tic; 
fprintf('CG for C_k: ')

tol=params.C.tol;
maxiter=params.C.maxiter;

L=Lambda2(F,G,Phi,abs(d)>0);


vec= @(x) x(:).'  ;
col= @(x) x(:);
ResC= @(x) reshape(x,size(C));
Resd= @(x) reshape(x,size(d));
ResB= @(x) reshape(x,size(B));

b=[col(d(:));col(sqrt(beta/2)*(B+Z./beta))];
init=col(C);

[x,flag,relres,iter,resvec]=lsqr(@afun,b,tol,maxiter,[],[],init); %add initial guess
    function y= afun(x,transp_flag)
        if strcmp(transp_flag,'notransp')
            y1=col(L*(ResC(x)));
            y2=sqrt(beta/2).*col(ResC(x));
            y=[y1;y2];
        else
            x1=x(1:numel(d));
            x2=x(numel(d)+1:end);
            y1=col(L'*(Resd(x1)));
            y2=sqrt(beta/2).*col(ResB(x2));
            y=[y1+y2];
        end
    end



Ck=ResC(x);

if params.visualization
figure(998);subplot(224);
plot((log10(abs(resvec)./norm(b(:)))),'r*-'); 
xlabel('iterations'); ylabel('10log of relative residual');
drawnow; end

t=toc; 
fprintf('t: %4.2f seconds',t)
fprintf('| relres %d | iters: %i | \n',relres,iter)

end
