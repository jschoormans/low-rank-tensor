function res = mtimes(a,b)
% res = mtimes(FT, x)
%
% reshape to 4D tensor
% ifft /fft 
% reshape back to 1-unfolded tensor

if a.adjoint 
res=reshape(b,a.tensorsize);
res=fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(res,1),2),[],1),[],2),1),2); %zero filled recon
res=res.*sqrt(a.size);
res=reshape(res,a.unfoldedsize);
else
res=reshape(b,a.tensorsize);
res=fftshift(fftshift(fft(fft(ifftshift(ifftshift(res,1),2),[],1),[],2),1),2); %zero filled recon
res=res./sqrt(a.size);
res=reshape(res,a.unfoldedsize);
end