% recon_13_08_VFA_fully sampled

folder=('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_08_13')
MR=MRecon([folder,'\lr_13082017_1643351_27_2_wip_sc23-vfa-t2prep_iV4.raw'])
%%
MR=MRecon([folder,'\lr_13082017_1643351_27_2_wip_sc23-vfa-t2prep_iV4.raw'])
MR.Parameter.Recon.ArrayCompression='no'; 
MR.Parameter.Recon.ACNrVirtualChannels=8;
MR.Parameter.Parameter2Read.chan=[34;35;36;37;44;45]
MR.ReadData;
MR.RandomPhaseCorrection;
MR.RemoveOversampling;
MR.PDACorrection;
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;
MR.SortData;
MR.GridData;
MR.RingingFilter;
MR.ZeroFill;

MR.Parameter.Recon.RemovePOversampling='No'
MR.RemoveOversampling
%%
%%
ifftc = @(x,n) ifftshift(ifft(fftshift(x,n),[],n),n);

sl=128;

ksp=MR.Data{1}; 
figure; plot(squeeze(abs(ksp(1,:,100,1,1,1))))
size(ksp);
ksp=ifftc(ksp,1); 
size(ksp)
ksp=ksp(:,65:65+127,:,:,:);
size(ksp)

che=create_checkerboard([1,1,size(ksp,3)]);
ksp=bsxfun(@times,ksp,che);
ksp=ksp(sl,:,:,:,:,:);
size(ksp)

im=ifftc(ifftc(ksp,2),3); 
size(im)

im=bart('rss 8',im);
size(im)
im=squeeze(im); 

figure(1); imshow(abs(im(:,:,1,1,1)),[])
figure(2); immontage4D(abs(im))
figure(3); immontage4D(angle(im),[-pi, pi])

%%
x=40; y=40; 
figure(4); 
plot(squeeze(abs(im(x,y,:,:,:)))); title(['plot over 5 dynamics of x=',num2str(x),' y=',num2str(y)]);

%% should fit T2 map 


tic
for xpix=1:128
    xpix
    parfor ypix=1:128
        
        xrange = [1:5].';
        yrange = squeeze(abs(im(xpix,ypix,:)));
        f = fit(xrange,yrange,'exp1');
        T2(xpix,ypix)=f.b;
        ampl(xpix,ypix)=f.a; 
    end
end
toc
%
figure(10); 

subplot(311); imshow(abs(im(:,:,1)),[]);
subplot(312); imshow(abs(T2),[]); 
colormap(jet)

%%
save('scan_27_sl_128','T2','im')


