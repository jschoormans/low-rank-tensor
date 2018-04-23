function N = gpuShift(M, k,varargin)
% shift 2D array on GPU
if nargin<3
    doubleoption=0;
else
    doubleoption=varargin{1};
end

if ~doubleoption
    N = zeros(size(M), 'single', 'gpuArray');
else
    N = zeros(size(M), 'double', 'gpuArray');
end

if k(1) >= 0  && k(2) >= 0
    N(k(1)+1:end, k(2)+1:end) = M(1:end-k(1), 1:end-k(2));
elseif k(1) >= 0  && k(2) < 0
    N(k(1)+1:end, 1:end+k(2)) = M(1:end-k(1), 1-k(2):end);
elseif k(1) < 0  && k(2) < 0
    N(1:end+k(1), 1:end+k(2)) = M(1-k(1):end, 1-k(2):end);
elseif k(1) < 0  && k(2) >= 0
    N(1:end+k(1), k(2)+1:end) = M(1-k(1):end, 1:end-k(2));
else
    error('wrong shift values');
end
end