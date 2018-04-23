function  res = Lambda2_hadam(F,G,Phi,mask,params)
% only for 2D

res.kspace_size = params.sizes.kspace;

res.F=F;
res.G=G;
res.Phi=Phi;
res.adjoint = 0;
res.mask=mask;
res = class(res,'Lambda2_hadam');

