
tol=1e-6;
maxiter=200;

L=Lambda(F,C,Phi,abs(du_1)>0);
a2PP=(alpha/2)*Psi'*Psi;

b=(L'*du_1 + (alpha/2)*Psi'*(A+Y./alpha));
b=b(:);
% try to reshape operator so we have proper matrix, vector calculations

X= @(G) (L'*(L*G)) +a2PP*G;
Res= @(x) reshape(x,[numel(G),1]);
ResA= @(x) reshape(x,size(G));
Aop = @(G) Res(X(ResA(G))); 

for i =1:40
    Aop(rand(size(G)));
end