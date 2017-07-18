function immontage4D(IM,varargin)
 %input: image tensor (4D)
 if nargin==2
     scale=varargin{1};
 else
     scale=[];
 end
 
Q=[];
for ii=1:size(IM,3)
    J=[];
    for jj=1:size(IM,4);
        J=[J,abs(IM(:,:,ii,jj))];
    end
    Q=[Q;J];
end

imshow(abs(Q),scale); axis off; 
end