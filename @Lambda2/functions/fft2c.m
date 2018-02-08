
function res = fft2c(x)
fctr = size(x,1)*size(x,2);
res = zeros(size(x));

% size_x = size(x);
% x = reshape(x, size_x(1), size_x(2), []);
% 
% for ii=1:size(x,3)
%         res(:,:,ii) = sqrt(fctr)*fftshift(ifft2(ifftshift(x(:,:,ii))));
% endx 
% res = (1/sqrt(fctr))*fftshift(fftshift(fft(fft(ifftshift(ifftshift(x,1),2),[],1),[],2),1),2);
res = (1/sqrt(fctr))*fftshift(fftshift(fft2(ifftshift(ifftshift(x,1),2)),1),2); %twice as fast (on GPU at least!) 

% res = reshape(res, size_x);

end