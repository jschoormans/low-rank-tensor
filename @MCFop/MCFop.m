function  res = MCFop(imsize,sens)
% only for 2D
res.sens=sens;
res.adjoint = 0;
res.imsize=imsize; 
res = class(res,'MCFop');

