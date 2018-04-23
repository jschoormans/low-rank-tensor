function [Kx2,Ky2]=findSharedKpoints_param2(mask,subspacedim)
%input: mask (kx,ky,subspacedim2)
% goal: find points which are sampled in all 1,dim2

Cdim2=squeeze(sum(mask(:,:,subspacedim,:),4))==size(mask,4);
assert(sum(Cdim2(:))>0,'no shared k-space points found for parameter dim 2!')

[Kx2,Ky2]=find(Cdim2);
end