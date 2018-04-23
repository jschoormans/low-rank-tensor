
clear

T2vals=sort([300,200,100,70,50,30]);



nsims=1e4; 

T2s=30+abs(randn(1,nsims)).^(1.3)*250;
figure(2); hist(abs(T2s),50); 
title('Monte Carlo T2 Distribution')

for i=1:nsims
    S(i,:)=rand*exp(-T2vals./T2s(i));
end
    
   figure(1); clf;hold on; 
[U,S,V]=svd(S,'econ')
for i=1:6
    plot(i+V(:,i)) 
end

title('SVD of T2 exponential decays')
S/sum(S(:))