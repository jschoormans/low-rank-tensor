function  res = Lambda(F,C,Phi,mask)
% only for 2D
res.F=F;
res.C=C;
res.Phi=Phi;
res.adjoint = 0;
res.mask=mask;
res = class(res,'Lambda');

