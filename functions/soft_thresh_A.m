function Ak=soft_thresh_A(G,Y,alpha,lambda,Psi)
% applies elementwise soft-thresholding operator
% as described in eqn B1
disp('soft-thresholding matrix A...')

T=Psi*G-Y./alpha;
% l2 norm of T
l2T=sqrt(sum(abs(T).^2,2));
thr=(l2T>lambda/alpha);

disp('---')
disp(['lambda/alpha: ',num2str(lambda/alpha)])
disp(['thresholding ',num2str(100*(sum(~thr(:))./numel(thr))),' percent'])
disp('---')

st=thr.*((l2T-(lambda/alpha))./(l2T+eps)); %soft thresholding
Ak=repmat(st,[1 size(T,2)]).*T;
end
