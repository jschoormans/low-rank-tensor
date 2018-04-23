function  res = TV_GPU(dim1,dim2,option,doubleoption)
% 2D Total variation operator that works on the GPU 

res.dim1=dim1;
res.dim2=dim2;
res.adjoint=0; %1= fft based , 2=circshift
% fft based is ~3 times faster, but is circular 


res.option=option;   
% pre-build convolution kernels;
zy=zeros(dim1,1); zy(1,1)=1; zy(2,1)=-1;
zx=zeros(1,dim2); zx(1,1)=1;  zx(1,2)=-1;
% zx=gpuArray(zx); zy=gpuArray(zy); 
zFx=fftn(zx); zFy=fftn(zy);
res.zFx=zFx; res.zFy=zFy; 

res.doubleoption=doubleoption; %1=double, 0=single 

res = class(res,'TV_GPU');

