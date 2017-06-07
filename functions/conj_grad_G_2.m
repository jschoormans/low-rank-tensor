function Gk=conj_grad_G_2(G,C,A,Y,alpha,Psi,du,Phi,F)
% new attempt
% minimizes the function
% argminG ||d - Fu G C Phi ||_2^2  + <Y,A-Psi G> + (alpha/2) ||A - Psi G ||_F ^2
% argminG ||d - Fu G C Phi ||_2^2  + <Y,A> -<Y,Psi G> + (alpha/2) ||A - Psi G ||_F ^2
% TO DO: CHECK GRADIENTS
% UNCOMMENT GRADIENT OF INNER PRODUCT 
% TUNE PARAMS
% MORE TESTING
% OPTIMIZE SPEED WITH PRECALCULATING FFTS AND ADDING SMARTER MASKING IN OBJ
% FUNTCTION 

% TEMP: 
Phi=Phi';

%params:
betals=0.6;
t0=1 ;
ls_alpha=1e-2;
maxlsiter=60;


iter=0;
s=0;
grad=calc_grad(G,C,Phi,du,F,Y,Psi,alpha); 

while(1)
    
    [f0,obj_l2,obj_inner_product,obj_F]= calc_objective(G,C,A,Y,alpha,Psi,du,Phi);

    % 1 calculate steepest direction
    gradk=calc_grad(G,C,Phi,du,F,Y,Psi,alpha);
    
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
        [f1]  =   calc_objective(G+t*sk,C,A,Y,alpha,Psi,du,Phi);
        lsiter=lsiter+1;
    end
    % update the position
    Gk=G+t*sk;
    
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
    G=Gk;
    iter=iter+1;
    if lsiter>=maxlsiter | (norm(s(:)) < 1e-8) 
        disp('CG convergenve reached...')
        break
    end
    
end
return


    function [obj,obj_l2,obj_inner_product,obj_F] = calc_objective(G,C,A,Y,alpha,Psi,du,Phi)
        % argminG ||d - Fu G C Phi ||_2^2  + <Y,A-Psi G> + (alpha/2) ||A - Psi G ||_F ^2
        du=reshape(du,[128^2 100]); %resize to 2D matrix TEMPORARY
        
        obj_l2_inner= du - (F*(G*C*Phi)) ;
        obj_l2_inner=obj_l2_inner.*(abs(du)>0); % make this more efficientS
        obj_l2=obj_l2_inner(:)'*obj_l2_inner(:);
        
        obj_inner_product=sum(sum(Y.*(A-Psi*G)));
        
        obj_F=sqrt(sum(sum(abs(A-Psi*G)).^2));
        
        obj=obj_l2+obj_inner_product+(alpha/2).*obj_F;
    end

    function gradk=calc_grad(G,C,Phi,du,F,Y,Psi,alpha)
        gl2=grad_l2(G,C,Phi,du,F);
        gip=grad_inner_product(Y,Psi);
        gf=grad_F(alpha,A,G,Psi);
        %         gradk=gl2+gip+gf; %temp: remove gip
        gradk=gl2+gf;

    end


    function grad=grad_l2(G,C,Phi,du,F)
        du=reshape(du,[128^2 100]); %resize to 2D matrix TEMPORARY 
%         grad = 2*Phi'*C'*(F'*(F*(G*C*Phi) -du)); %du is a vector --> shd be a tensor'/ 
        grad = 2* (F'*(du * Phi'*C' + F*(G*C*Phi*Phi'*C'))); %new attempt...
    end

    function grad = grad_inner_product(Y,Psi)
        % <y, Psi * G>
        grad= -2*(Y'*Psi); % not sure at all..
    end

    function grad= grad_F(alpha,A,G,Psi)
        grad= alpha*(2*Psi'*A+G); 
    end

end
