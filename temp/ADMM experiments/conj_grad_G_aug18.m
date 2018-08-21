function Gk=conj_grad_G_aug18(G,C,A,Y,alpha,Psi,d,Phi,F,params)
%implement reduced problem size

tic
samplingmask=abs(d)>0;
tol=params.G.tol;
maxiter=params.G.maxiter;


%input data

vec= @(x) x(:).'  ;
col= @(x) x(:);

ResG= @(x) reshape(x,size(G));
Resd= @(x) reshape(x,size(d));
ResA= @(x) reshape(x,size(A));


    function y= afun2(x,transp_flag)
        if strcmp(transp_flag,'notransp')
            y1=col(samplingmask.*(F*(ResG(x)*C*Phi)));
            y2=sqrt(alpha/2).*col(Psi*(ResG(x)));
            y=[y1;y2];
        else
            x1=x(1:numel(d));
            x2=x(numel(d)+1:end);
            y1=col((F'*(Resd(x1))*Phi'*C'));
            y2=sqrt(alpha/2).*col(Psi'*ResA(x2));
            y=[y1+y2];
        end
    end


CGoption=params.G.CGoption

if CGoption==1
    
    b=[col(d(:));col(sqrt(alpha/2)*(A+Y./alpha))];
    init=col(G);
    
    [x,flag,relres,iter,resvec]=lsqr(@afun2,b,tol,maxiter,[],[],init); %add initial guess
    Gk=ResG(x);
    
    if params.visualization
        figure(998);subplot(223);
        plot((log10(abs(resvec)./norm(b(:)))),'r*-');
        xlabel('iterations'); ylabel('10log of relative residual');
        drawnow; end
    
    
elseif CGoption==2
    
    %use lsqr algo - but with ATA as input
    ATA2=@(x) col(F'*(samplingmask.*(F*ResG(x)*C*Phi))*Phi'*C'+(alpha/2).*(Psi'*(Psi*ResG(x))));

    rhs = col((F'*(Resd(col(d(:))))*Phi'*C')) + sqrt(alpha/2).*col(Psi'*ResA(col(sqrt(alpha/2)*(A+Y./alpha))));
    init=col(G);
    
    [x,flag,relres,iter,resvec]=cgs(ATA2,rhs,tol,maxiter,[],[],init); %add initial guess
    Gk=ResG(x);
    
    if params.visualization
        figure(998);subplot(223);
        plot((log10(abs(resvec)./norm(rhs(:)))),'r*-');
        xlabel('iterations'); ylabel('10log of relative residual');
        drawnow; end
    
    
else
    
    ATA=@(x) (F'*(samplingmask.*(F*ResG(x)*C*Phi))*Phi'*C')+(alpha/2).*(Psi'*(Psi*ResG(x)));
    x=G;
    rhs = col((F'*(Resd(col(d(:))))*Phi'*C')) + sqrt(alpha/2).*col(Psi'*ResA(col(sqrt(alpha/2)*(A+Y./alpha))));
    
    r=ResG(rhs)-ATA(x);
    
    p=r;
    k=1;
    r2=r(:).'*r(:);
    
    clear CGresidual
    while k<60 %&& r2 > 1e-4
        k
        alphaCG=r2/((conj(vec(p(:)))*col(ATA(p))));
        x=x+alphaCG.*p;
        
        rold=r; r2old=r2;
        r=r-alphaCG.*ATA(p);
        
        r2=r(:).'*r(:);
        
        beta=r2/(r2old);
        p=r+beta.*p;
        k=k+1;
        CGresidual(k)=r2;
        
    end
    figure(10); clf;
    if exist('CGresidual')
        semilogy(abs(CGresidual)); title('CG step (update x)'); drawnow;
    end
    
    
    Gk=ResG(x);
    
end



%%


t=toc;
fprintf('t: %4.2f seconds',t)
% fprintf('| relres %d | iters: %i | \n',relres,iter)
end