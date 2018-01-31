
res1=200; res2=150; 

Psi1=opConvolve(res1,res2,[-1 1],[0 0],'cyclic');
Psi2=opConvolve(res1,res2,[-1 1]',[0 0],'cyclic');
Psi=[Psi1;Psi2];


TVG=TV_GPU(200,150);

R=phantom(200);
R=R(:,1:150);

RG=gpuArray(R);

tic;
for ii=1:10;
temp1 = TVG*RG(:);
temp3= TVG'*temp1(:);

end
toc


tic;
for ii=1:10;
temp2 = Psi*R(:);
temp4 = Psi'*temp2(:);

end
toc


figure(1); subplot(211); imshow(abs(cat(2,reshape(gather(temp1),[200 300]),reshape(temp2,[200 300]))),[0 0.1])
figure(1); subplot(212); imshow(abs(cat(2,R,reshape(gather(temp3),[200 150]),reshape(temp4,[200 150]))),[0 4])




