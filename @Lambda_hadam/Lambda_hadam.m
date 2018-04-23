function  res = Lambda_hadam(F,C,Phi,mask,params)
% only for 2D
% lambda operator: input G, output: Omega(F*G*C*Phi)
% 6/7/2017: included pre-calculated matrix products CxPhi and Phi*xC* 

% Added by Bobby 31-01-2018
res.kspace_size = params.sizes.kspace;
%
res.F=F;
% res.C=C;
% res.Phi=Phi;
res.PhiTCT=Phi'*C';     % precalc to save time 
res.CPhi=C*Phi;         % pre calc to save time 
res.adjoint = 0;
res.mask=mask;
res = class(res,'Lambda_hadam');

