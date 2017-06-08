function Bk=soft_thresh_B(C,Z,mu,beta)
mub=mu/beta;
CZb=C-(Z./beta); %precalculate to speed up 

thr=abs(CZb)>=mub; %threshold

disp('---')
disp(['mu/beta: ',num2str(mub)])
disp(['thresholding ',num2str(100*(sum(~thr(:))./numel(thr))),' percent'])

Bk=thr.*(abs(CZb)-mub).*((CZb)./(abs(CZb)+eps));
end