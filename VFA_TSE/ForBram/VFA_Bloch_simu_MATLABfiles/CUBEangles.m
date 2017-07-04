function [Angles] = CUBEangles(T1,T2,TE,Star,etl)

%calculate flip angles for CUBE from F(-1/1) and Z(-1/1)
%output: Array with size [1,etl] with flip angles
%input: T1 (ms), T2 (ms), TE (ms), Starget (between 0-1) and ETL

%matrices for extended phases of  F (Mxy) and Z  (Mz)
%number of Fs and Zs
nrmagn = 3+(etl-1)*4;
F=zeros(nrmagn,etl);
Z=zeros(nrmagn,etl);
s=zeros(1,etl);

%index i of F(i) -> i+round(nrmagn/2)
ishift = round(nrmagn/2);

%start magnetization
F(0+ishift,1)=1;
Z(0+ishift,1)=0;

%half TE for precession
tau=TE/2;

for j=1:etl
    
    G=F; %temp variables for loop
    Y=Z;
    
    for f=2:nrmagn %precession after echo
        F(f,j)=G(f-1,j)*exp(-tau/T2);
        Z(f,j)=Y(f,j)*exp(-tau/T1);
    end
    
    if j==1
        s(j)=0.5*(1+Star);
    else
        s(j)=0.5*(s(j-1)+Star);
    end
    
    Angles(j)=CUBEangle(Z(1+ishift,j),F(-1+ishift,j),F(1+ishift,j),T2,s(j),tau);
    
    G=F; %temp variables for loop
    Y=Z;
    
    for f=1:nrmagn %rotation F(i) = F(f), F(-i) = F(nrmagn-f+1)
        [F(f,j),Z(f,j)] = CUBEnutation(G(nrmagn-f+1,j),G(f,j),Y(f,j),Angles(j));
    end
    
    G=F; %temp variables for loop
    Y=Z;
    
    for f=2:nrmagn %precession before echo
        F(f,j)=G(f-1,j)*exp(-tau/T2);
        Z(f,j)=Y(f,j)*exp(-tau/T1);
    end
    
    %pass magn on to next echo
    if j<etl
        F(:,j+1) = F(:,j);
        Z(:,j+1) = Z(:,j);
    end
       
end

end
 
function [angle] = CUBEangle(Zp,Fm,Fp,T2,St,tau) %Z(+1),F(-1),F(1)
    sqtrm = Zp^2-(Fp-St*exp(tau/T2))*(Fm-St*exp(tau/T2));
    if sqtrm>=0
        angle=2*atan( (Zp + sqrt(sqtrm)) / (Fp - St*exp(tau/T2)));
    else
        angle = pi;
    end
end

function [Fkn, Zkn] = CUBEnutation(Fm,Fp,Zp,a) %F(-i) and F(i)!
    [NM] = [0.5*(1+cos(a)) 0.5*(1-cos(a)) sin(a); -0.5*sin(a) 0.5*sin(a) cos(a)] * [Fp;Fm;Zp];
    Fkn = NM(1);
    Zkn = NM(2);
end





