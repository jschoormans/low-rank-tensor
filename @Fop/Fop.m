function  res = Fop(tensorsize,unfoldedsize)
% only for 2D
res.adjoint = 0;
res.tensorsize=tensorsize;
res.unfoldedsize=unfoldedsize;
res.size=tensorsize(1)*tensorsize(2);
res = class(res,'Fop');

