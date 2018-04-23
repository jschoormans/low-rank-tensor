function [Kx1,Ky1]=findSharedKpoints_param1(mask,subspacedim)
%input: mask (kx,ky,subspacedim1)
% goal: find points which are sampled in all dim1,1

%temporary change Bobby 16-01-2018. Changed back 22-01-2018.
Cdim1=squeeze(sum(mask(:,:,:,subspacedim),3))==size(mask,3);
% Cdim1=squeeze(sum(mask(:,:,:),3))==size(mask,3);

assert(sum(Cdim1(:))>0,'no shared k-space points found for parameter dim 1!')

[Kx1,Ky1]=find(Cdim1);
end