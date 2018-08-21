
% try low-rank plus sparse simulation...


N=64
T2vals=[100,100,100,100,100,100,100,100]
T1vals=[500,700,900,10,700,880,200,105]
T2prep=[10,20,30,40,50,100,200,600]
TI=[100]
I=T2_T1_phantom(N,T2vals,T1vals,T2prep, TI,1,0); % low rank
I=I+double(rand(size(I))>0.995); % add sparse info


imagine(squeeze(abs(I)))

%% define functions and operators

nt=size(I,4)
rt=@(x) reshape(x,[N^2, nt])
r2t=@(x) reshape(x,[N,N, nt])
rproj=@(x) reshape(x,[N,N*nt])
vec = @(x) x(:);

shrinkage=@(x,lambda) (x./(abs(x)+eps)).*max(abs(x)-lambda,0);

F=opDFT2(N,N,1);
F=opBlockDiag(F,F,F,F,F,F,F,F)

allpatterns=[]
for ii=1:nt
    pattern=((rand(N,N)>0.75))==1;
    pattern(25:40,25:40)=ones(16);
    allpatterns=[allpatterns(:);pattern(:)];
end
E=opExcise(F,squeeze(allpatterns)==1,'rows')

% construct k-space
k=E*I(:);  % kspace

%% undersample k-space
M=E'*k;
imagine(r2t(abs(M)))

%% TV operator..

T=opConvolve(N^2,nt,[1 -1],[0 0],'regular'); %unsure about cyclic of course...

Tx=opConvolve(N,N,[1 -1].',[0 0],'regular'); %unsure about cyclic of course...
Ty=opConvolve(N,N,[1 -1],[0 0],'regular'); %unsure about cyclic of course...
T1=[Tx;Ty]
T=opBlockDiag(T1,T1,T1,T1,T1,T1,T1,T1)

TM=T*M(:);


invT=pinv(T)
%
lambda_S=0.1
lambda_L=1
niter=100
visualize=1

M0=E'*k;

M=M0;
S=zeros(size(M));
L=zeros(size(M));

l2norm=zeros(niter,1);

clear nuclear_norm l1norm
for iter=1:niter
    if visualize==1
        figure(1);
        imshow(cat(1,abs(rproj(M)),abs(rproj(M0))),[0 1]); title(['iter=',num2str(iter)]); drawnow;
    end
    fprintf('iter %f\n',iter)
    fprintf('update L...\n')
    Lold=L;
    [LSVT,SS]=SVT(rt(M-S),lambda_L);
    L=vec(LSVT);
    
    fprintf('update S...\n')
    S=invT*(shrinkage(T*(M-Lold),lambda_S));
    
    fprintf('Data consistency...\n')
    M=L+S-E'*(E*(L+S)-k);
    
    l2norm(iter)=norm(E*(L+S)-k,2);
    l1norm(iter)=lambda_S*norm(T*S,1);
    nuclear_norm(iter)=lambda_L*sum(abs(SS(:)));
    
    fprintf('l2-norm=%4.2f...\n',l2norm(iter))
    fprintf('l1-norm=%4.2f...\n',l1norm(iter))
    fprintf('nuclear-norm=%4.2f...\n',nuclear_norm(iter))
    
    if visualize==1
        
        figure(2);
        imshow(cat(1,abs(rproj(L)),abs(rproj(S))),[0 1]);
        title(['L(upper), S(lower), iter=',num2str(iter)]); drawnow;
        figure(3);
        plot(l1norm); hold on; plot([1:niter],l2norm,'-+'); plot(nuclear_norm); hold off
        legend('l1','l2','nuclear norm')
    end
end

