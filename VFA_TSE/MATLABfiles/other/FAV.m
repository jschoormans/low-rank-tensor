tr=14;
fa=[1:0.1:90];
fat=fa*pi/180;
At=1000;
T1 = 700
R1t = 1/T1;


yd=At*(1-exp(-tr*R1t)).*sin(fat) ./ (1-exp(-tr*R1t).*cos(fat));
yd=yd+0.09*randn(size(yd));
scatter(fa,yd,'.k')
