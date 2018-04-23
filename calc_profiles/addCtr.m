function mask= addCtr(mask, ctrsize)
%input: ctrsize =-1; - do not add center

assert(ndims(mask)==2); 
[ky,kz]=size(mask);
kyc=floor((1+ky)/2);
kzc=floor((1+kz)/2);

if ctrsize~=-1
mask(kyc-ctrsize:kyc+ctrsize,kzc-ctrsize:kzc+ctrsize)=ones(2*ctrsize+1,2*ctrsize+1); 
else 
   % do nothing  
end
    
end