T1 = 1200;
T2 = 38;
TE = 5.2;
Star = 0.1;
etl = 67;
k=1;
j=1;
lastangle=0;

while lastangle < pi && k < 30
    
    Star = Star + 0.01;
    angles = CUBEangles(T1,T2,TE,Star,etl);
	lastangle = angles(end);
    k=k+1;
end
    
Star = Star-0.01;
angles = CUBEangles(T1,T2,TE,Star,etl);
lastangle = angles(end)

while lastangle < pi && j < 30
    
    Star = Star + 0.001;
    angles = CUBEangles(T1,T2,TE,Star,etl);
	lastangle = angles(end);
    j=j+1;
end

Star = Star-0.001;
angles = CUBEangles(T1,T2,TE,Star,etl);
lastangle = angles(end);

while lastangle < pi && j < 30
    
    Star = Star + 0.0001;
    angles = CUBEangles(T1,T2,TE,Star,etl);
	lastangle = angles(end);
    j=j+1;
end

% T2 = T2;
Star = Star-0.0001;
angles = CUBEangles(T1,T2,TE,Star,etl);
lastangle = angles(end)
Star


