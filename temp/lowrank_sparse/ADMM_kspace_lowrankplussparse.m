% m = 64*64;
% n = 50;

% r = 4;
% 
% L = randn(m,r) * randn(r,n);    % low rank
% S = sprandn(m,n,0.05);          % sparse
% S(S ~= 0) = 20*binornd(1,0.5,nnz(S),1)-10;
% V = 0.01*randn(m,n);            % noise
%%

Nx=64
T2vals=[100,500,600,10,30,50,80,100]
T1vals=[500,700,900,10,700,880,200,105]
T2prep=[10,20,30,40,50,100,200,600]
TI=[100]
L=T2_T1_phantom(Nx,T2vals,T1vals,T2prep, TI,1,0); % low rank
S=double(rand(size(I))>0.995); % add sparse info
V=zeros(size(S));
A = S + L + V;

% imagine(squeeze(abs(A)))

%% define functions and operators

nt=size(I,4)
rt=@(x) reshape(x,[Nx^2, nt])
r2t=@(x) reshape(x,[Nx,Nx, nt])
rproj=@(x) reshape(x,[Nx,Nx*nt])
vec = @(x) x(:);



%%
N = 3; %??

A=rt(I);
m=size(A,1)
n=size(A,2) 
F=opDFT2(Nx,Nx,0);
% E=opBlockDiag(F,F,F,F,F,F,F,F)


allpatterns=[]
for ii=1:nt
    pattern=((rand(Nx,Nx)>0.02))==1;
    pattern(25:40,25:40)=ones(16);
    figure(3); imshow(pattern);
    F1{ii}=opExcise(F,~pattern,'rows')

end
% allpatterns=ones(size(allpatterns)); allpatterns(1)=0;
E=opBlockDiag(F1{1},F1{2},F1{3},F1{4},F1{5},F1{6},F1{7},F1{8})

% construct k-space
kspace=E*I(:);  % kspace

Tx=opConvolve(Nx,Nx,[1 -1].',[0 0],'regular'); %unsure about cyclic of course...
Ty=opConvolve(Nx,Nx,[1 -1],[0 0],'regular'); %unsure about cyclic of course...
T1=[Tx;Ty]
T=opBlockDiag(T1,T1,T1,T1,T1,T1,T1,T1)
invT=pinv(T); 

figure(1); imshow(A,[])
%%
g2_max = norm(A(:),inf);
g3_max = norm((T*A(:)));
g2 = 2e-2*g2_max;  % WHAT ARE THESE???
g3 = 0.01*g3_max;


MAX_ITER = 100;
ABSTOL   = 1e-4;
RELTOL   = 1e-1;

tic;

lambda = 1;
rho = 1/lambda;


m=size(kspace,1)

X_1 = zeros(m,n);
X_2 = zeros(m,n);
X_3 = zeros(m,n);
z   = zeros(m,N*n);
U   = zeros(m,n);

fprintf('\n%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
    'r norm', 'eps pri', 's norm', 'eps dual', 'objective');

for k = 1:MAX_ITER

    B = avg(X_1, X_2, X_3) - rt(kspace)./N + U;

    % x-update
    X_1 = (1/(1+lambda))*(X_1 - B);
%     X_2 = rt(T'*(prox_l1(T*(X_2(:) - B(:)), lambda*g2)));
    X_2 = rt(invT*(prox_l1(T*(X_2(:) - B(:)), lambda*g2)));

    X_3 = prox_matrix(X_3 - B, lambda*g3, @prox_l1);

    % (for termination checks only)
    x = [X_1 X_2 X_3];
    zold = z;
    z = x + repmat(-avg(X_1, X_2, X_3) + rt(kspace)./N, 1, N);

    % u-update
    U = B;

    % diagnostics, reporting, termination checks
    h.objval(k)   = objective(X_1, g2, X_2, g3, X_3);
    h.r_norm(k)   = norm(x - z,'fro');
    h.s_norm(k)   = norm(-rho*(z - zold),'fro');
    h.eps_pri(k)  = sqrt(m*n*N)*ABSTOL + RELTOL*max(norm(x,'fro'), norm(-z,'fro'));
    h.eps_dual(k) = sqrt(m*n*N)*ABSTOL + RELTOL*sqrt(N)*norm(rho*U,'fro');

    if k == 1 || mod(k,10) == 0
        fprintf('%4d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
            h.r_norm(k), h.eps_pri(k), h.s_norm(k), h.eps_dual(k), h.objval(k));
    end

    if h.r_norm(k) < h.eps_pri(k) && h.s_norm(k) < h.eps_dual(k)
         break;
    end

end

h.admm_toc = toc;
h.admm_iter = k;
h.X1_admm = X_1;
h.X2_admm = X_2;
h.X3_admm = X_3;


figure(2); 
imshow(cat(1,rproj(abs(X_1)),rproj(abs(X_2)),rproj(abs(X_3))),[])

fprintf('\nADMM (vs true):\n');
fprintf('|V| = %.2f;  |X_1| = %.2f\n', norm(rt(V), 'fro'), norm(X_1,'fro'));
fprintf('nnz(S) = %d; nnz(X_2) = %d\n', nnz(rt(S)), nnz(X_2));
fprintf('rank(L) = %d; rank(X_3) = %d\n', rank(rt(L)), rank(X_3));
