function [Starg,lastangle] = CUBElastangle_v3(T1,T2,dTE,Stargini,ETLe,maxangle)

maxangle=maxangle*(pi/180);
Starg = Stargini;
k=1;
j=1;
lastangle=0;

while lastangle < maxangle && k < 30
    
    Starg = Starg + 0.01;
    angles = CUBEangles(T1,T2,dTE,Starg,ETLe);
	lastangle = angles(end);
    k=k+1;
end
    
Starg = Starg-0.01;
angles = CUBEangles(T1,T2,dTE,Starg,ETLe);
lastangle = angles(end);

while lastangle < maxangle && j < 30
    
    Starg = Starg + 0.001;
    angles = CUBEangles(T1,T2,dTE,Starg,ETLe);
	lastangle = angles(end);
    j=j+1;
end

Starg = Starg-0.001;
angles = CUBEangles(T1,T2,dTE,Starg,ETLe);
lastangle = angles(end);

while lastangle < maxangle && j < 30
    
    Starg = Starg + 0.0001;
    angles = CUBEangles(T1,T2,dTE,Starg,ETLe);
	lastangle = angles(end);
    j=j+1;
end

% T2 = T2;
Starg = Starg-0.0001;
angles = CUBEangles(T1,T2,dTE,Starg,ETLe);
lastangle = angles(end);
