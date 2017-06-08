function res = mtimes(a,b)
% res = mtimes(FT, x)
F=a.F;
C=a.C;
Phi=a.Phi;

if a.adjoint %L'(d) ===> G
    res=F'*(b)*Phi'*C'; % not totally sure about this operator.
else %L*G= F G C Phi
    res = F*b*C*Phi;
    res=res.*a.mask;

end