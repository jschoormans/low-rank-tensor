function  res = MCFop(imsize,sens)
% multi-coil DFT operator 
res.sens=sens;
res.ncoils=size(sens,3);
res.adjoint = 0;
res.imsize=imsize; 
res = class(res,'MCFop');

