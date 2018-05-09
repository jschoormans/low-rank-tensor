function res = mtimes(a,bb)


if (a.adjoint)
    if ndims(bb)==2
        reshapeflag=1;
        sbb=size(bb); 
        bb= reshape(bb,[size(bb,1),1,a.dim1,a.dim2]); %for precon CG (matrix input)
    else
        reshapeflag=0;
    end

    
    
    for ii=1:a.dim1
        for jj=1:a.dim2
            res(:,:,:,:,ii,jj)=a.O{ii,jj}'*(gather(bb(:,:,ii,jj)).*sqrt(a.w(:,ii,jj)));
        end
    end
     if reshapeflag==1
        res=reshape(res,[size(res,1)*size(res,2)*size(res,3)*size(res,4),size(res,5)*size(res,6)]) ;
     end
    
else
    
    
    if ndims(bb)==2
        reshapeflag=1;
        sbb=size(bb);
        bb= reshape(bb,[size(bb,1),1,1,1,a.dim1,a.dim2]); %for precon CG (matrix input)
    else
        reshapeflag=0;
    end
    
    
    for ii=1:a.dim1
        for jj=1:a.dim2
            res(:,:,ii,jj)=a.O{ii,jj}*gather(bb(:,:,:,:,ii,jj));
        end
    end
    
    if reshapeflag==1
        res=reshape(res,[size(res,1)*size(res,2),size(res,3)*size(res,4)]) ;
    end
    
    
end