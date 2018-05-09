function  res = Lambda(F,C,Phi,mask)
% only for 2D
% lambda operator: input G, output: Omega(F*G*C*Phi)
% 6/7/2017: included pre-calculated matrix products CxPhi and Phi*xC* 

res.F=F;
% res.C=C;
% res.Phi=Phi;
res.PhiTCT=Phi'*pinv(C);     % precalc to save time  %changed 9-5-2018
res.CPhi=C*Phi;         % pre calc to save time 
res.adjoint = 0;
res.mask=mask;
res = class(res,'Lambda');

