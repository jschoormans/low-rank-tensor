function wmatf=getRadWeightsGA(ku)
nx=size(ku,1);
nspokes=size(ku,2);
for ii=1:size(ku,3); %for different bins/frames/etc...


angles=(360*(atan2(imag(ku(1,:,ii)),real(ku(1,:,ii)))+pi))./(2*pi); %angle of each spoke
[anglesup,anglesindex]=sort(mod(angles,180).*(pi)/180);
delta_anglesup=anglesup-circshift(anglesup,1,2);
delta_anglesup=mod(delta_anglesup,pi);
R=1; %temp
wr=abs(2*R*tan(delta_anglesup/2))./(2*pi); %division for normalization
wmean=(wr+circshift(wr,-1,2))./2;

wmax=nx*(pi/2)./nspokes;

wvec=abs(2*(0:nx-1)/nx-1);

wmat(:,:,ii)=wvec'*wmean*wmax*nspokes;
wmat(nx/2+1,:,ii)=ones(1,nspokes).*(pi/(4*nspokes));
wmatf(:,anglesindex,ii)=wmat(:,:,ii);
end
end