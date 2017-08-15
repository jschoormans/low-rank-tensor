function Ak=soft_thresh_A(G,Y,alpha,lambda,Psi,operatorsize)
% applies elementwise soft-thresholding operator
% as described in eqn B1

T=Psi*G-Y./alpha;
% l2 norm of T
l2T=sqrt(sum(abs(T).^2,2));
thr=(l2T>lambda/alpha);

fprintf('soft-thresholding A: ')
fprintf('lambda/alpha: %1.2e  |',(lambda/alpha))
fprintf('thresholding %1.2f  percent \n',(100*(sum(~thr(:))./numel(thr))))

st=thr.*((l2T-(lambda/alpha))./(l2T+eps)); %soft thresholding
Ak=repmat(st,[1 size(T,2)]).*T;

title_text = sprintf('s.t. A: %d data points thresholded (%1.2f%%).',sum(~thr(:)), (100*(sum(~thr(:))./numel(thr))));
figure(39); clf; spy(reshape(thr,[operatorsize])); title(title_text); 
end
