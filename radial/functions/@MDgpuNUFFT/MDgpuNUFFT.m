function [res] = MDgpuNUFFT(k,w,osf,wg,sw,imageDim,sens,varargin)
% for now, keep w and sens indepdent of dimensions... to do: change this 
kDims=size(k); 

if ndims(k)>2
    res.dim1=kDims(3);
else 
    res.dim1=1; 
end
if ndims(k)>3
res.dim2=kDims(4);
else 
    res.dim2=1; 
end

for ii=1:res.dim1
    for jj=1:res.dim2
        %make all operators...
        O{ii,jj}= gpuNUFFT(k(:,:,ii,jj),w(:,ii,jj),osf,wg,sw,imageDim,sens,true,true,true);
    end
end

res.O=O; 
res.w=w; 
res.adjoint = 0;

res = class(res,'MDgpuNUFFT');
