%% @2Dt (eg diffusion /DCE whatever)000000000000000000

res=64;
A=phantom(res);
B=exp(-linspace(0,2,res));
C=repmat(B,[res 1 res]); C=permute(C,[1 3 2]);
R=rand(res,res)>0.5; R=repmat(R,[1 1 res]);
D=repmat(A,[1 1 res]);

E=C+C.*D.*R;


%%
 %tensor ranks
a=10
b=a;
c=2;

[F1,F2,F3,F4]=tucker(E,[a,b,c]);

disp('fraction of original information needed:')
res^3./(res*a+res*b+res*c+a*b*c)
disp('pct of orignal signal explained:' ); F3
