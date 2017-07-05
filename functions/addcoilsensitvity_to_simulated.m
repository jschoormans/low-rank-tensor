%%  add simulated coil sensitivity 
function [Iout,S_normalized]=addcoilsensitvity_to_simulated(I,ncoils)
[nx ny ndim1,ndim2]=size(I);
assert(nx==ny)

S=bart(['phantom -S',num2str(ncoils),' -x',num2str(nx)]);
S=squeeze(S); 
S_norm=sqrt(sum(abs(S).^2,3));
S_normalized=bsxfun(@rdivide,S,S_norm);

Iout= bsxfun(@times,permute(I,[1 2 5 3 4]),conj(S_normalized));


% Iout=permute(Iout,[1 2 5 3 4]);
end