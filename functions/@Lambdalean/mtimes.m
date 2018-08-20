function res = mtimes(a,b)
% res = mtimes(FT, x)
F=a.F;

if a.adjoint %L'(d) ===> G
    %     res = F'*(b*a.PhiTCT);
    %     res = (F'*b)*a.PhiTCT;
%     
%     parfor n=1:150
%     Q=a.Phi(:,n)*a.Phi(:,n).'
%     end
    res = (F'*(b))*a.PhiTCT;
    
else %L*G= F G C Phi
    
    
    res = F*(b*a.CPhi);
    res = res.*a.mask;
    
end