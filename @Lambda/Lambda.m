function  res = Lambda(F,C,Phi)
% only for 2D
res.F=F;
res.C=C;
res.Phi=Phi;
res.adjoint = 0;
res = class(res,'Lambda');

