% recon 28 july
clear all; close all; clc
% cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_07_28')  
cd('/home/qzhang/lood_storage/divi/Ima/parrec/Jasper/LRT/Low_Rank_2017_09_21');
%
% MR=MRecon('lr_28072017_1118286_15_2_wip_vfa-t2prep_tenr18dynsV4.raw'); % wrong orientation/resolution/t2prepvals
% MR1=MRecon('lr_28072017_1242319_23_2_wip_vfa-t2prep_tenr18dynsV4.raw')
% MR1=MRecon('lr_28072017_1112088_14_2_wip_vfa-t2prep_tenr18dynsV4.raw')
% MR1.Perform
%%
MR=MRecon('lr_21092017_1843595_6_2_wipsc23vfat2preptenr18dynsV4.raw')  %OIL
MR=MRecon('lr_21092017_2002407_19_2_wipsc23vfat2preptenr18dynsV4.raw') %BRAIN
DTI=false;

MR.Parameter.Parameter2Read.typ = 1;


% load data
disp('readdata')
MR.ReadData;
MR.RandomPhaseCorrection;
disp('corrections...')
MR.RemoveOversampling;
MR.PDACorrection; %???
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;

disp('sortdata')

MR.Parameter.Labels.Index.aver=MR.Parameter.Labels.Index.rf;
MR.Parameter.Recon.ImmediateAveraging='No'
MR.SortData;

K=MR.Data;
K=ifftshift(ifft(K,[],1),1); size(K)
K=K(40,:,:,:,:,:,:,:,:,:,:,:);
% remove stupid checkerboard pattern
che=create_checkerboard([1,size(K,2),size(K,3)]);
K=bsxfun(@times,K,che);
a = sum(K(:,:,:,:,1,:,1,1,1,1),6)./sum(K(:,:,:,:,1,:,1,1,1,1)~=0,6);
sens=bart('ecalib -S -m1',a);
% sens=ones(size(kspace));
kspace=squeeze(K);


if DTI
par_dim1 = 9; par_dim2 = 6;
du = zeros(size(kspace,1),size(kspace,2),size(kspace,3),par_dim1,par_dim2);
for ii = 1:par_dim1
    for jj = 1:par_dim2
        du(:,:,:,ii,jj) = kspace(:,:,:,((jj-1) * par_dim1 +ii));
    end
end
kspace=du; 
end
%%

if ndims(kspace)==4 % if one chan
   kspace=permute(kspace,[1 2 5 3 4]); 
end

kspace=kspace(2:end-1,:,:,:,:);
a = sum(sum(kspace,4),5)./sum(sum(kspace~=0,4),5);
sens=bart('ecalib -S -m1',permute(a,[4 1 2 3]));
sens=sens+1e-7; % no zero vals in sense maps...

sens=sens;
close all;
params=params_init();
params.Lg=4;
params.inspectLg=false
params.L3=4;
params.L4=4;
params.sparsity_transform='TV';
params.Imref=[];
params.x=17;
params.y=62;
params.lambda=2e-3
% sens(sens==0)=1e-2;

params.increase_penalty_parameters=false
params.G.precon=false;
params.G.maxiter=20
params.niter=3;

coil=2
P_recon=LRT_recon(kspace(:,:,:,:,:),squeeze(sens(:,:,:,:)),params);
% visualize recon
figure(1000); immontage4D(squeeze(abs(P_recon)),[]);
%%
x1=20
x2=33
x3=40
y1=36
y2=40
y3=33

figure(200);clf;subplot(221); hold on 
imshow(squeeze(abs(P_recon(:,:,1,5,1))),[]);
plot(y1,x1,'r+')
plot(y2,x2,'b+')
plot(y3,x3,'g+'); hold off
title('first echo, lowest TE');

subplot(223); hold on 
imshow(squeeze(abs(P_recon(:,:,1,1,20))),[]);
plot(y1,x1,'r+')
plot(y2,x2,'b+')
plot(y3,x3,'g+'); hold off
title('last echo, highest TE');
subplot(222)
hold on
plot(squeeze(abs(P_recon(x1,y1,1,:,3))),'r')
plot(squeeze(abs(P_recon(x2,y2,1,:,3))),'b')
plot(squeeze(abs(P_recon(x3,y3,1,:,3))),'g')
hold off
title('T2prep dimension')

subplot(224)
hold on
plot(squeeze(abs(P_recon(x1,y1,1,1,:))),'r')
plot(squeeze(abs(P_recon(x2,y2,1,1,:))),'b')
plot(squeeze(abs(P_recon(x3,y3,1,1,:))),'g')
hold off
title('echo dimension')
