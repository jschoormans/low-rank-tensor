function Bk=soft_thresh_B(C,Z,mu,beta)
mub=mu/beta;
CZb=C-(Z./beta); %precalculate to speed up 

thr=abs(CZb)>=mub; %threshold
Bk=thr.*(abs(CZb)-mub).*((CZb)./(abs(CZb)+eps));
end