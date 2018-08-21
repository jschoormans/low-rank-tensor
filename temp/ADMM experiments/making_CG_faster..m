

% option email 
% not correct and not fast... 

tic 
FGC=F*(G*C); 
parfor ii=1:150
    t(:,:,ii)=samplingmask(:,ii).*(FGC*(Phi(:,ii)*Phi(:,ii)'));
end

a=sum(t,3);
a=F'*(a*C')
toc

% option before 
tic
b=(F'*(samplingmask.*(F*ResG(G)*C*Phi))*Phi'*C');
toc