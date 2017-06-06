% simulation of a DCE - inverison recovery measurement 
res=64
nt=100; 
nTI=40; 

A=phantom(res);

A=repmat(A,[1 1 nt nTI]);

B=1-2*exp(-linspace(0,2,nTI));
B=permute(B,[1 3 4 2]); % set b dimenion too 4(nti)
B=repmat(B,[res res nt 1]);

% C: time dimensinos
C=rand(1,1,nt,1);
C=repmat(C,[res res 1 nTI]);

% make image;

D=A.*B.*C;

size(D)
size(C)
size(B)
size(A)

%% TUCKER DEMCOPOSITION

E=reshape(D,[res.^2, nt nTI]); % group dims

E=E+randn(size(E)).*0.001;
a=10
b=10
c=10

[F1,F2,F3,F4]=tucker(E,[a b c]);

F3
%%


disp('fraction of original information needed:')
(res^2*nt*nTI)./(res^2*a+nt*b+nTI*c+a*b*c)
disp('pct of orignal signal explained:' ); F3
