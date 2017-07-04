%%  add simulated coil sensitivity 
function [Iout,S_normalized]=addcoilsensitvity_to_simulated(I,ncoils)
[nx ny ndim1,ndim2]=size(I);
assert(nx==ny)
x1=[linspace(-1,1,nx)];
[meshgrid1,meshgrid2]=meshgrid(x1,x1);


muvec=linspace(-0.7,0.7,ncoils);
muvec2=0.7*(-1).^(1:ncoils);

% Iout=[];
for i =1:ncoils
    Si=mvnpdf([meshgrid1(:),meshgrid2(:)],[muvec(i) muvec2(i)],[.5 0.2 ; 0.2 0.5]); % 2D Gaussian 
    Si=reshape(Si,[nx ny]);
%     Iout=cat(5,Iout,I.*repmat(Si,[1 1 ndim1 ndim2]));
    S(:,:,i)=Si;
end
S_norm=sqrt(sum(abs(S).^2,3));
S_normalized=bsxfun(@rdivide,S,S_norm);

Iout= bsxfun(@times,permute(I,[1 2 5 3 4]),S_normalized);


% Iout=permute(Iout,[1 2 5 3 4]);
end