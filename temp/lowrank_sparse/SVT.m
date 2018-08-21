function [xt,Ss]=SVT(x,lambda) 

shrinkage=@(x) (x./(abs(x)+eps)).*max(abs(x)-lambda,0);
[U,S,V] = svd(x,'econ');
Ss=shrinkage(S);
xt=U*Ss*V;

end