function Gk=conj_grad_G_2(G,C,Ak,Bk,Y,Z,alpha,Psi,du,Phi,P0_1)
% new attempt 
% minimizes the function 
% argminG ||d - Fu G C Phi ||_2^2  + <Y,A-Psi G> + (alpha/2) ||A - Psi G ||_F ^2
% argminG ||d - Fu G C Phi ||_2^2  + <Y,A> -<Y,Psi G> + (alpha/2) ||A - Psi G ||_F ^2 


obj = calc_objective(G,C,Ak,Y,alpha,Psi,du,Phi)



end

function obj = calc_objective(G,C,A,Y,alpha,Psi,du,Phi)
% argminG ||d - Fu G C Phi ||_2^2  + <Y,A-Psi G> + (alpha/2) ||A - Psi G ||_F ^2 


obj_l2_inner= du - F*G*C*Phi ; % is Phi okay here?? % to do ; F operator 
obj_l2=obj_l2_inner.'*obj_l2_inner;

obj_inner_product=sum(sum(Y.*(A-Psi*G)));

obj_F=sqrt(sum(sum(abs(A-Psi*G)).^2));

obj=obj_l2+obj_inner_product+(alpha/2).*obj_F;
end


function grad=grad_l2(G,C,Phi,du,F)
grad = 2*F'*(F*G*C*Phi -du);
end

function grad = grad_inner_product(Y,Psi)
% <y, Psi * G>

grad= Psi'*Y;                %% ?? really unsure about this gradient?? 
end

function grad= grad_F(alpha,A,G,Psi)

grad_F= (alpha)* Psi' * (A-Psi * G);
end


