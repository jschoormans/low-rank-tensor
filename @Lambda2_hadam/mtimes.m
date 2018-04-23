function res = mtimes(a,b)
F=a.F;
G=a.G;
Phi=a.Phi;
kspace_size = a.kspace_size;
psi_hadam = 0.5*[[1,1,1,1];[1,1,-1,-1];[1,-1,1,-1];[1,-1,-1,1]];
reshape_1 = [kspace_size(1)*kspace_size(2)*kspace_size(4),kspace_size(5)];
reshape_2 = [kspace_size(1)*kspace_size(2),kspace_size(4)*kspace_size(5)];

if a.adjoint %L'(d) ===> C
    res = (F'*(b));
    % Now preperation for hadamard transform and reshape and transpose back.
    res=reshape(((psi_hadam*((reshape(res,reshape_1)).')).'),reshape_2);
    %
    res=G'*((res)*Phi'); % not totally sure about this operator. (stond hier al, van Jasper/Kerry?)
    
    % It is assumed that the inverse of psi_hadam is itself (which is right
    % I think).
else %L*C= F G C Phi
    res = (G*b*Phi);
    % Now preperation for hadamard transform and reshape and transpose back.   
    res = reshape(((psi_hadam*((reshape(res,reshape_1)).')).'),reshape_2); % Same as adjoint case above.
    %
    res = F*(res);
    res=res.*a.mask;
end