function PDF = generatePDF(m,ctrsize,nr_points)
%make probability density function 
%input m: (matrix of mask size)
% ctrsize, size of fully sampled center
% nr points: total number of points


m1=linspace(-1,1,size(m,1)); 
m2=linspace(-1,1,size(m,2)); 

coords=sqrt(bsxfun(@plus,(m1.').^2,m2.^2));
coords=(sqrt(2)-(coords))./sqrt(2); 
coords=coords.^2; 

PDF=zeros(size(m)); 
PDF=addCtr(PDF,ctrsize);

% calculate total area of PDF 
nr_centerpoints=sum(PDF(:)); 
nr_points_outside=nr_points-nr_centerpoints;

PDFoutside=coords./((sum(sum(coords(~PDF))))/nr_points_outside); %
PDF=max(PDF,PDFoutside);
PDF=1-PDF;
end