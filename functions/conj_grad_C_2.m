function Gk=conj_grad_G_2(G,C,Y,B,Z,beta,du,Phi,F)
% new attempt
% minimizes the function
% argminG ||d - Fu G C Phi ||_2^2  + <Z,B-Z> - (beta/2) ||B -C ||_F ^2


%params:
betals=0.6;
t0=1 ;
ls_alpha=1e-2;
maxlsiter=60;

iter=0;
s=0;
grad=calc_grad(G,C,Phi,du,F,Y,Psi,alpha); 

while(1)
    
    [f0,obj_l2,obj_inner_product,obj_F]= calc_objective(F,G,C,Phi,B,Z,du,beta);

    % 1 calculate steepest direction
    gradk=calc_grad();
    
    % 2 compute beta
    beta=(gradk(:).'*gradk(:))/(grad(:).'*grad(:)+eps); %Fletcher-Reeves;
    
    % 3 update conjugate direction
    sk=-grad+ beta*s;
    
    % 4 perform a line search :
    t=t0;
    lsiter=0;
    f1=1e99; %only to start...
    while (f1 > f0 - ls_alpha*t*abs(gradk(:)'*sk(:)))^2 & (lsiter<maxlsiter)
        t = t * betals;
        [f1]  =   calc_objective(F,G,C+t*sk,Phi,B,Z,du,beta);
        lsiter=lsiter+1;
    end
    % update the position
    Ck=C+t*sk;
    
    % print some iteration comments
    disp(['iter: ',num2str(iter),'| lsiter: ',num2str(lsiter), '| obj:',num2str(f0)])

    %update parameters for next iteration;
    if lsiter > 2
		t0 = t0 * betals;
	end 
	if lsiter<1
		t0 = t0 / betals;
	end

    grad=gradk;
    s=sk;
    C=Ck;
    iter=iter+1;
    if lsiter>=maxlsiter | (norm(s(:)) < 1e-8) 
        disp('CG convergenve reached...')
        break
    end
    
end
return


    function [obj,obj_l2,obj_inner_product,obj_F] = calc_objective(F,G,C,Phi,B,Z,du,beta)
        du=reshape(du,[128^2 100]); %resize to 2D matrix TEMPORARY
        
        obj_l2_inner= du - (F*(G*C*Phi)) ;
        obj_l2_inner=obj_l2_inner.*(abs(du)>0); % make this more efficientS
        obj_l2=obj_l2_inner(:)'*obj_l2_inner(:);
        
        obj_inner_product=trace(Z'*(B-Z));
        
        obj_F=norm((B-C),'fro')^2 ;
        
        obj=obj_l2+obj_inner_product-(beta/2).*obj_F;
    end

    function gradk=calc_grad()
        gl2=grad_l2();
        gip=grad_inner_product();
        gf=grad_F();
        gradk=gl2+gip+gf;

    end


    function grad=grad_l2()
    end

    function grad = grad_inner_product()
    end

    function grad= grad_F()
    end

end
