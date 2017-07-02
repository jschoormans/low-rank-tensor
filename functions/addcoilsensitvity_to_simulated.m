%%  add coil sensitivity information 
function [Iout,S]=addcoilsensitvity_to_simulated(I,ncoils)
[nx ny ndim1,ndim2]=size(I)
assert(nx==ny)
x1=[linspace(-1,1,nx)]
[meshgrid1,meshgrid2]=meshgrid(x1,x1)


muvec=linspace(-0.6,0.6,ncoils)
Iout=[];
for i =1:ncoils
    Si=mvnpdf([meshgrid1(:),meshgrid2(:)],[muvec(i) 0],[.5 0.2 ; 0.2 0.5]);
    Si=reshape(Si,[nx ny]);
    Iout=cat(5,Iout,I.*repmat(Si,[1 1 ndim1 ndim2]));
    S(:,:,i)=Si;
end


Iout=permute(Iout,[1 2 5 3 4]);
end