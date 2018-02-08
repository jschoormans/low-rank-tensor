
function res = fft2c_fftmod(x)
% same time as fftshift essentially, but only works for even kspace input.
% could be faster if the multplications are combined with other steps in
% the operators....

fctr = size(x,1)*size(x,2);
res = zeros(size(x));

if mod(size(x,1),4)==0
    fftmod1=((-1).^[0:size(x,1)-1]);
elseif mod(size(x,1),2)==0
    fftmod1=((-1).^[1:size(x,1)]).*1i;
else
    error('N should be even!')
end

if mod(size(x,2),4)==0
    fftmod2=((-1).^[0:size(x,2)-1]);
elseif mod(size(x,1),2)==0
    fftmod2=((-1).^[1:size(x,2)]).*1i;
else
    error('N should be even!')
end

res=bsxfun(@mtimes,bsxfun(@mtimes,x,fftmod1.'),fftmod2);
res = fft2(res);
res=bsxfun(@mtimes,bsxfun(@mtimes,res,fftmod1.'),fftmod2);
res=(1/sqrt(fctr))*res;

% res = (1/sqrt(fctr))*fftshift(fftshift(fft2(ifftshift(ifftshift(x,1),2)),1),2); %twice as fast (on GPU at least!) 

end