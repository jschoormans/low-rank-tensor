function res = mtimes(a,b)
% res = mtimes(TV_GPU, x)

% TV*x - where x is a vector of size [res1,res2]
if a.adjoint %TV'*x
    % reshape source vector
    res=reshape(b,[a.dim1,a.dim2,2]);
    res1=bsxfun(@ldivide,fftn(res(:,:,1)),a.zFx);
    res2=bsxfun(@ldivide,fftn(res(:,:,2)),a.zFy);
    res=ifftn(res1)+ifftn(res2);
    res=res(:); 
    
    
else  %TV*x 
    % input vector b, output: TVx,TVy
    % reshape source vector 
    res=reshape(b,[a.dim1,a.dim2]);
    res=fftn(res);
    
    %TVx 
    res1=bsxfun(@mtimes,a.zFx,res);
    res1=ifftn(res1);
    
    %TVy
    res2=bsxfun(@mtimes,a.zFy,res);
    res2=ifftn(res2);

    res=[res1(:);res2(:)];
    
    
end