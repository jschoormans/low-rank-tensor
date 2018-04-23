function res = mtimes(a,b)
% res = mtimes(FT, x)
F=a.F;
kspace_size = a.kspace_size;
psi_hadam = 0.5*[[1,1,1,1];[1,1,-1,-1];[1,-1,1,-1];[1,-1,-1,1]];
reshape_1 = [kspace_size(1)*kspace_size(2)*kspace_size(4),kspace_size(5)];
reshape_2 = [kspace_size(1)*kspace_size(2),kspace_size(4)*kspace_size(5)];

if a.adjoint %L'(d) ===> G
% res = (F'*(b))*a.PhiTCT; %original
res = (F'*(b)); % stap 1: outpout P (b=kspace)
% Now preperation for hadamard transform and reshape and transpose back.
res=reshape(((psi_hadam*((reshape(res,reshape_1)).')).'),reshape_2);
%
res=res*a.PhiTCT; % stap 2: output G    (P=GCPhi; G=P Phi^T C^T )

% It is assumed that the inverse of psi_hadam is itself (which is right I
% think).
else %L*G= F G C Phi
    res = b*a.CPhi;
% Now preperation for hadamard transform and reshape and transpose back.   
    res = reshape(((psi_hadam*((reshape(res,reshape_1)).')).'),reshape_2); % Same as adjoint case above.
%
    res = F*res;    
    res = res.*a.mask;
end