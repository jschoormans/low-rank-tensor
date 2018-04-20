function res = mtimes(a,b)
% res = mtimes(TV_GPU, x)
% note: tried diff(X) --> just as fast as circshift
% todo: vectorizing (removing loop) might make this faster
ncols=size(b,2); 
% further optimization by removing the for loop for pagefun???

% TV*x - where x is a vector of size [res1,res2]

if a.adjoint %TV'*x
    
    if a.option>3 %dont use loop options
        
        if a.option==4
            rest=reshape(b,[a.dim1,a.dim2,2,ncols]);
            res1=pagefun(@plus,...
                pagefun(@min,rest(:,:,1,:),gpuShift(rest(:,:,1,:),[0 -1],a.doubleoption)),...
                pagefun(@min,rest(:,:,2,:),gpuShift(rest(:,:,2,:),[-1 0],a.doubleoption)));
            res=reshape(res1,a.dim1*a.dim2,ncols);
        end
        
    else % use loop options(deprecated - only use for CPU!!!)
        for coliter=1:ncols
            
            if a.option==1; %%FFT BASED
                % reshape source vector
                rest=reshape(b(:,coliter),[a.dim1,a.dim2,2]);
                res1=bsxfun(@rdivide,fftn(rest(:,:,1)),(a.zFx+eps));
                res2=bsxfun(@rdivide,fftn(rest(:,:,2)),(a.zFy+eps));
                rest=ifftn(res1)+ifftn(res2);
                res(:,coliter)=rest(:);
            elseif a.option==2 % normal circshift
                rest=reshape(b(:,coliter),[a.dim1,a.dim2,2]);
                res1=rest(:,:,1)-circshift(rest(:,:,1),-1,2)+rest(:,:,2)-circshift(rest(:,:,2),-1,1);
                res(:,coliter)=res1(:);
            elseif a.option==3 %option 3GPU circshift
                rest=reshape(b(:,coliter),[a.dim1,a.dim2,2]);
                res1=rest(:,:,1)-gpuShift(rest(:,:,1),[0 -1],a.doubleoption)+rest(:,:,2)-gpuShift(rest(:,:,2),[-1 0],a.doubleoption);
                res(:,coliter)=res1(:);
                
            end
            
            
            
        end
    end
    
    
else  %TV*x 
    % input vector b, output: TVx,TVy
    % reshape source vector
    
        if a.option>3 %dont use loop options
        
        if a.option==4
            rest=reshape(b,[a.dim1,a.dim2,ncols]);
            
            res(:,:,1,:)=pagefun(@min,rest,gpuShift(rest,[0 1],a.doubleoption));
            res(:,:,2,:)=pagefun(@min,rest,gpuShift(rest,[1 0],a.doubleoption));
            res=reshape(res,[a.dim1*a.dim2*2,ncols]);
        end
        else
    
    
    
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
        elseif a.option==2 % circshift 
            res1=rest-circshift(rest,1,2);
            res2=rest-circshift(rest,1,1);
            res(:,coliter)=[res1(:);res2(:)];
            
        elseif a.option==3 %option 3: gpuShift 
            res1=rest-gpuShift(rest,[0 1],a.doubleoption);
            res2=rest-gpuShift(rest,[1 0],a.doubleoption);
            res(:,coliter)=[res1(:);res2(:)];

        end
        
    end
        end  
    end

end
