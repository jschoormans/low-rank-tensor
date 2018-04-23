%monte carlo making probability density function!

nr_MCpoints=nr_points-nr_centerpoints;
pdf1=bsxfun(@times,ones(size(m)),permute(linspace(-1/sqrt(2),1/sqrt(2),size(m,1)).^2,[2 1]));
pdf2=bsxfun(@times,ones(size(m)),permute(linspace(-1/sqrt(2),1/sqrt(2),size(m,2)).^2,[1 2]));
pdf=1-pdf1-pdf2;
pdf=addCtr(pdf,ctrsize);
scaling=nr_MCpoints/(sum(pdf(pdf~=1)));
pdf=pdf.*scaling;
surf(pdf)
pdf=addCtr(pdf,ctrsize);
% sum(pdf(:))
% nr_MCpoints
% nr_points