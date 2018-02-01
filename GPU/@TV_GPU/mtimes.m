function res = mtimes(a,b)
% res = mtimes(TV_GPU, x)
% note: tried diff(X) --> just as fast as circshift
% todo: vectorizing (removing loop) might make this faster
ncols=size(b,2); 


% TV*x - where x is a vector of size [res1,res2]
if a.adjoint %TV'*x
    
    for coliter=1:ncols
        
        if a.option==1; %%FFT BASED
            % reshape source vector
            rest=reshape(b(:,coliter),[a.dim1,a.dim2,2]);
            res1=bsxfun(@rdivide,fftn(rest(:,:,1)),(a.zFx+eps));
            res2=bsxfun(@rdivide,fftn(rest(:,:,2)),(a.zFy+eps));
            rest=ifftn(res1)+ifftn(res2);
            res(:,coliter)=rest(:);
        else
            rest=reshape(b(:,coliter),[a.dim1,a.dim2,2]);
            res1=rest(:,:,1)-circshift(rest(:,:,1),-1,2)+rest(:,:,2)-circshift(rest(:,:,2),-1,1);
            res(:,coliter)=res1(:);
        end
        
        
        
    end
    
    
else  %TV*x 
    % input vector b, output: TVx,TVy
    % reshape source vector
    for coliter=1:ncols
        
        rest=reshape(b(:,coliter),[a.dim1,a.dim2]);
        
        if a.option==1; % FFT BASED
            rest=fftn(rest);
            
            %TVx
            res1=bsxfun(@mtimes,a.zFx,rest);
            res1=ifftn(res1);
            
            %TVy
            res2=bsxfun(@mtimes,a.zFy,rest);
            res2=ifftn(res2);
            
            res(:,coliter)=[res1(:);res2(:)];
        else
            res1=rest-circshift(rest,1,2);
            res2=rest-circshift(rest,1,1);
            res(:,coliter)=[res1(:);res2(:)];
            
        end
        
        
        
    end

end