function  res = TV_GPU(dim1,dim2)
% 2D Total variation operator that works on the GPU 

res.dim1=dim1;
res.dim2=dim2;
res.adjoint=0; 

% pre-build convolution kernels;
zy=zeros(dim1,1); zy(1,1)=1; zy(2,1)=-1;
zx=zeros(1,dim2); zx(1,1)=1;  zx(1,2)=-1;
% zx=gpuArray(zx); zy=gpuArray(zy); 
zFx=fftn(zx); zFy=fftn(zy);
res.zFx=zFx; res.zFy=zFy; 


res = class(res,'TV_GPU');

