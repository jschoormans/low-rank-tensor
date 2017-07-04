function [Mss,Mz] = CUBEmagn(T1,T2,TE,Angles,etl,StMagn)

%Calculates extended phase graph magnetizations

%output: 
% Mss = Transversal magnetization phase graph. Middle row is phase 0. 
% Mz = Longitudinal magnetization phase graph. Middle row is phase 0.
%input: 
% T1 (ms), T2 (ms), TE (ms), Flip angles (rad), etl, Start Magn [0,1]

%matrices for F and Z to calculate their values prior to RF 

%number of Fs and Zs
nrmagn = 3+(etl-1)*4;
F=zeros(nrmagn,etl);
Z=zeros(nrmagn,etl);

%index i of F(i) -> i+round(nrmagn/2)
ishift = round(nrmagn/2);

%start magnetization
F(0+ishift,1)=StMagn;
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

Mss=F;
Mz=Z;

end

function [Fkn, Zkn] = CUBEnutation(Fm,Fp,Zp,a) %F(-i) and F(i)!
    [NM] = [0.5*(1+cos(a)) 0.5*(1-cos(a)) sin(a); -0.5*sin(a) 0.5*sin(a) cos(a)] * [Fp;Fm;Zp];
    Fkn = NM(1);
    Zkn = NM(2);
end





