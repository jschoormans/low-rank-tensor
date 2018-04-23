function arraymasked = multiplymask(array,mask)
% applies a mask to a given array. The mask is a matrix and multiplication
% is done along the first two dimensions.

arrayunfolded = reshape(array,[size(array,1), numel(array)/size(array,1)]);
repmask = repmat(mask,1,numel(array)/(size(array,1)*size(array,2)));
arraymaskedunfolded = arrayunfolded.*repmask;
arraymasked = reshape(arraymaskedunfolded,size(array));