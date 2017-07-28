function sortArray(MRD,ky,kz,chan,param1)
ky=MR.Parameter.Labels.Index.ky;
kz=MR.Parameter.Labels.Index.kz;
chan=MR.Parameter.Labels.Index.chan;
param1=MR.Parameter.Labels.Index.dyn;
MRD=MR.Data;

nkx=size(MR.Data,1)

minky=min(ky)
maxky=max(ky)
nky=length(minky:maxky)
kyunique=unique(ky)

minkz=min(kz)
maxkz=max(kz)
nkz=length(minky:maxky)
kzunique=unique(kz)

minchan=min(chan)
maxchan=max(chan)
chanunique=unique(chan)
nchans=length(chanunique)

p1unique=unique(param1)
minp1=min(p1unique)
maxp1=max(p1unique)
np1=length(p1unique);

D=zeros(nkx,nky,nkz,nchans,np1);

kyindex=ky-minky+1;
kzindex=kz-minkz+1;
chanindex= chan-1; % should be able to deal with eg chan=2,3,5,6
np1index= param1+1; % TO DO: CHANGE: 

tic
for ii=1:size(MRD,2)
    if mod(ii,1e3)==0;waitbar(ii/size(MRD,2));end
    D(:,kyindex(ii),kzindex(ii),chanindex(ii),np1index(ii))=MR.Data(:,ii);
end
toc

end