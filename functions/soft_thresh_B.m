function Bk=soft_thresh_B(C,Z,mu,beta)
mub=mu/beta;
CZb=C-(Z./beta); %precalculate to speed up 

thr=abs(CZb)>=mub; %threshold


fprintf('soft-thresholding B: ')
fprintf('mu/beta: %1.2e  |',(mub))
fprintf('thresholding %1.2f  percent \n',(100*(sum(~thr(:))./numel(thr))))


Bk=thr.*(abs(CZb)-mub).*((CZb)./(abs(CZb)+eps));

title_text = sprintf('s.t. C: %d data points thresholded (%1.2f%%).',sum(~thr(:)), (100*(sum(~thr(:))./numel(thr))));
figure(49);clf; spy(thr); title(title_text); 
end