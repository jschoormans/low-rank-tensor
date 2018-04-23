%
% size(K) =     1   249    66     2     3    50



%%
Kc1=squeeze(K(1,:,:, 1,1,:)); 
size(Kc1)
mask=squeeze(sum(abs(Kc1)>0,3));

mask=zeros(size(mask)); 

ctrsize=3; 
mask(125-floor(ctrsize/2):125+floor(ctrsize/2),33-floor(ctrsize/2):33+floor(ctrsize/2))=50*ones(ctrsize,ctrsize);
size(mask);

idx=find(mask==50);

Kc1_reshaped=reshape(Kc1,[size(Kc1,1)*size(Kc1,2),size(Kc1,3)]);

InputMatrix=Kc1_reshaped(idx,:); 
size(InputMatrix);

figure(2); imshow(InputMatrix,[]); colormap jet
%
clear S V D
[S V D]=svd(InputMatrix.','econ')

figure(3); plot(abs(S(:,3)),'r')