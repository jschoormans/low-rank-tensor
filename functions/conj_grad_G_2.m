function Gk=conj_grad_G_2(G,C,A,Y,alpha,Psi,du,Phi,F)
% new attempt
% minimizes the function
% argminG ||d - Fu G C Phi ||_2^2  + <Y,A-Psi G> + (alpha/2) ||A - Psi G ||_F ^2
% argminG ||d - Fu G C Phi ||_2^2  + <Y,A> -<Y,Psi G> + (alpha/2) ||A - Psi G ||_F ^2


%params:
beta=0.6;
t0=1 
ls_alpha=1;
maxlsiter=40;


iter=0
s=0;
grad=calc_grad(G,C,Phi,du,F,Y,Psi,alpha); 

while(1)
    
    f0 = calc_objective(G,C,A,Y,alpha,Psi,du,Phi);
    
    % 1 calculate steepest direction
    gradk=calc_grad(G,C,Phi,du,F,Y,Psi,alpha);
    
    % 2 compute beta
    beta=(gradk.'*gradk)/(grad.'*grad+eps); %Fletcher-Reeves;
    
    % 3 update conjugate direction
    sk=-grad+ beta*s;
    
    % 4 perform a line search :
    t=t0;
    lsiter=0;
    while (f1 > f0 - ls_alpha*t*abs(gradk(:)'*sk(:)))^2 & (lsiter<maxlsiter)
        lsiter = lsiter + 1;
        t = t * beta;
        [f1]  =   calc_objective(G,C,A,Y,alpha,Psi,du,Phi);
    end
    
    
    % update the position
    Gk=G+t*sk;
    
    % print some iteration comments
    disp(['iter: ',num2str(iter),' lsiter: ',num2str(lsiter), ' '])
    %
    
    %update parameters for next iteration;
    grad=gradk;
    s=sk;
    
    iter=iter+1;
    if iter>itermax % TO DO: add stopping criterion for convergence
        break
    end
    
end

    function obj = calc_objective(G,C,A,Y,alpha,Psi,du,Phi)
        % argminG ||d - Fu G C Phi ||_2^2  + <Y,A-Psi G> + (alpha/2) ||A - Psi G ||_F ^2
        
        obj_l2_inner= du - F*G*C*Phi ; % is Phi okay here?? % to do ; F operator
        obj_l2=obj_l2_inner.'*obj_l2_inner;
        
        obj_inner_product=sum(sum(Y.*(A-Psi*G)));
        
        obj_F=sqrt(sum(sum(abs(A-Psi*G)).^2));
        
        obj=obj_l2+obj_inner_product+(alpha/2).*obj_F;
    end

    function gradk=calc_grad(G,C,Phi,du,F,Y,Psi,alpha)
        gl2=grad_l2(G,C,Phi,du,F);
        gip=grad_inner_product(Y,Psi);
        gf=grad_F(alpha,A,G,Psi);
        gradk=gl2+gip+gf;
    end


    function grad=grad_l2(G,C,Phi,du,F)
        Phi=Phi.';
        du=reshape(du,[128^2 100]); %resize to 2D matrix TEMPORARY 
        grad = 2*Phi'*C'*(F'*(F*(G*C*Phi) -du)); %du is a vector --> shd be a tensor'/ 
    end

    function grad = grad_inner_product(Y,Psi)
        % <y, Psi * G>
        grad= Psi'*Y;                %% ?? really unsure about this gradient??
    end

    function grad= grad_F(alpha,A,G,Psi)
        grad_F= (alpha)* Psi' * (A-Psi * G);
    end

end
