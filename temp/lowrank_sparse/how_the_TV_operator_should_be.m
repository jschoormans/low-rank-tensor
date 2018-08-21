
T=opConvolve(N^2,nt,[-1 1],[0 0],'cyclic'); %unsure about cyclic of course...

lambda_S=0.1
S=pinv(T)*(shrinkage(T*(M0),lambda_S));
figure(1);
imshow(cat(1,abs(rproj(S)),abs(rproj(M0))),[0 1]); title(['iter=',num2str(iter)]); drawnow;
