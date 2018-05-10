function Gk=conj_grad_G_new(G,C,A,Y,alpha,Psi,d,Phi,F,params)
tic
samplingmask=abs(d)>0; 
tol=params.G.tol;
maxiter=params.G.maxiter;

L=Lambda(F,C,Phi,samplingmask);

%input data

vec= @(x) x(:).'  ;
col= @(x) x(:);

ResG= @(x) reshape(x,size(G));
Resd= @(x) reshape(x,size(d));
ResA= @(x) reshape(x,size(A));

    
    
b=[col(d(:));col(sqrt(alpha/2)*(A+Y./alpha))];
init=col(G);;

[x,flag,relres,iter,resvec]=lsqr(@afun2,b,tol,maxiter,[],[],init); %add initial guess
    function y= afun2(x,transp_flag)
        if strcmp(transp_flag,'notransp')
            y1=col(L*(ResG(x)));
            y2=sqrt(alpha/2).*col(Psi*(ResG(x)));
            y=[y1;y2];
        else
            x1=x(1:numel(d));
            x2=x(numel(d)+1:end);
            y1=col(L'*(Resd(x1)));
            y2=sqrt(alpha/2).*col(Psi'*ResA(x2));
            y=[y1+y2];
        end
    end
Gk=ResG(x);

if params.visualization
figure(998);subplot(223);
plot((log10(abs(resvec)./norm(b(:)))),'r*-'); 
xlabel('iterations'); ylabel('10log of relative residual');
drawnow; end

t=toc; 
fprintf('t: %4.2f seconds',t)
fprintf('| relres %d | iters: %i | \n',relres,iter)
end
