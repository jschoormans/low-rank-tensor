function res = mtimes(a,b)
% res = mtimes(FT, x)
%
% reshape to 5D tensor of size (imsize(1),imsize(2),nc,s4,s5);
% ifft /fft 
% reshape back to 1-unfolded tensor

if size(b,1)==a.imsize(1)&& size(b,2)==a.imsize(2)
    if a.adjoint %F'*B where B is a tensor of size (imsize(1),imsize(2),nc,s4,s5);
        res=fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(b,1),2),[],1),[],2),1),2); %zero filled recon
        res=res.*sqrt(a.imsize(1)*a.imsize(1));
        res=bsxfun(@mtimes,res,a.sens)
        
    else %F'*B where B is a tensor of size (imsize(1),imsize(2),nc,s3,s4);
        res=fftshift(fftshift(fft(fft(ifftshift(ifftshift(b,1),2),[],1),[],2),1),2); %zero filled recon
        res=res./sqrt(a.imsize(1)*a.imsize(2));
        res=sum(bsxfun(@mtimes,res,conj(a.sens)),3);

    end
    
    
elseif size(b,1)==a.imsize(1)*a.imsize(2) %reshape matrix before fft
    
    res=reshape(b,[a.imsize(1),a.imsize(2),numel(b)/(a.imsize(1)*a.imsize(2))]);
    
    if a.adjoint %F'*B where B is a tensor of size (imsize(1),imsize(2),s3,s4);
        res=fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(res,1),2),[],1),[],2),1),2); %zero filled recon
        res=res.*sqrt(a.imsize(1)*a.imsize(1));
        res=res.*sqrt(a.imsize(1)*a.imsize(1));
        res=bsxfun(@mtimes,res,a.sens)

    else %F'*B where B is a tensor of size (imsize(1),imsize(2),s3,s4);
        res=fftshift(fftshift(fft(fft(ifftshift(ifftshift(res,1),2),[],1),[],2),1),2); %zero filled recon
        res=res./sqrt(a.imsize(1)*a.imsize(2));
        res=sum(bsxfun(@mtimes,res,conj(a.sens)),3);

    end
    
    res=reshape(res,size(b));
else
    error('unsupported matrix size!')
end