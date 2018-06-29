clear; close all; clc
cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2018_01_31_VFA')
% recon T2-VFA copied from: recon_10_31_modfrom_5_9.m


%%

fn='lr_31012018_1733043_2_2_wip_brin-vfa-t2prep-lowrankV4.raw'

clear MR
MR=MRecon(fn); 

DTI=0;
% MR.Parameter.Parameter2Read.chan=[10; 11];
% MR.Parameter.Parameter2Read.ky=[-80:80].';

MR.Parameter.Labels.Index.aver=(MR.Parameter.Labels.Index.rf);
MR.Parameter.Parameter2Read.aver=[0:max(MR.Parameter.Labels.Index.aver)].'; 
MR.Parameter.Parameter2Read.typ = 1;
MR.Parameter.Recon.ImmediateAveraging='No';
% MR.Parameter.Parameter2Read.chan=[34;35;36;37;44;45]
% load data
disp('readdata')
MR.ReadData;
MR.RandomPhaseCorrection;
disp('corrections...')
MR.RemoveOversampling;
MR.PDACorrection; 
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;
MR.K2IM;

%% from here onwards we recon only one slice

slice_nr=100; 

%   ifft in readout direction + select slice 
%   K=fftshift(ifft(ifftshift(MR.Data,1),[],1),1);
MR2=MR.Copy;
MR2.Data=MR.Data(slice_nr,:,:,:,:,:);
K=sortArray(MR2);
if true
    disp('performing coil compresion...')
    Kcc=bart('cc -p5',permute(K,[1 2 3 4 7 8 9 10 5 6]));
else
    Kcc=K;
end

Kcc=permute(Kcc,[1 2 3 4 9 10 5 6 7 8]);
size(Kcc)
%%
kspace=Kcc(1,1:end,1:end,:,:,:,:,:,:,:,:,:);
imshow(squeeze(abs(kspace(1,:,:,1,1))),[])
size(kspace)

% remove stupid checkerboard pattern
che=create_checkerboard([1,size(kspace,2),size(kspace,3)]);
kspace=bsxfun(@times,kspace,che);
kspace=squeeze(kspace);
a = sum(sum(kspace(:,:,:,:,:),4),5)./(sum(sum(kspace(:,:,:,:,:)~=0,4),5)+eps);
% a2 = (sum(kspace(:,:,:,7,:),5))./((sum(kspace(:,:,:,7,:)~=0,5))+eps);

sens=bart('ecalib -S -r15 -m1',permute(a,[4 1 2 3]));
% sens=sens+1e-7; % no zero vals in sense maps...

figure(112)
subplot(311)
immontage4D(abs(sens),[])
subplot(312)
immontage4D(real(sens),[])
subplot(313)
immontage4D(angle(sens),[-pi pi])

%%
kspaceinput=kspace(:,:,:,:,3:end); %dummy pulses 

params=params_init();
params.L3=3;
params.L4=3;
params.subspacedim1=2;
params.subspacedim2=3; 
params.scaleksp=1; 

params.Lg=8;
params.inspectLg=false;
params.sparsity_transform='TV';


params.Imref=[];
params.x=15;
params.y=120;
params.mu=0.3e3;
params.lambda=6e-2;
params.alpha=2
params.beta=2
params.autolambda=0; 
params.automu=1; 
% sens(sens==0)=1e-2;

params.niter=10; 
params.increase_penalty_parameters=false; 
params.G.precon=true;
params.G.maxiter=50;



P_recon=LRT_recon_test(kspaceinput(:,:,1,:,:),squeeze(sens(:,:,:,1)),params);
%% visualize recon
figure(1000); immontage4D(squeeze(abs(P_recon)),[]);

%% make gif for all contrasts 
cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2018_01_31_VFA\recons_may18')
for jj=1:6
h=  figure(1001); clf;
set(gcf,'Color','White')
filename=['gif-contrast-GPU',num2str(jj),'.gif']
for i=1:48; 
    imshow(squeeze(abs(P_recon(:,:,:,jj,i))).',[],'InitialMagnification',100,'Border','Tight'); drawnow;
    
          frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 

          % Write to the GIF File 
      if i == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1); 
      end 

    
end

end
%% direct iFFT 

recon=bart('fft -i 3',kspace); 
recon=bart('rss 4',recon); 

figure(1001); immontage4D(squeeze(abs(recon)),[]);


%%
figure(1002); 
subplot(121); imshow(cat(2,squeeze(abs(P_recon(:,:,:,1,1))),squeeze(abs(P_recon(:,:,:,1,5)))),[]);
subplot(122); imshow(cat(2,squeeze(abs(recon(:,:,:,1,1))),squeeze(abs(recon(:,:,:,1,5)))),[])

%%
figure(1003); 
imshow(rot90(squeeze(abs(P_recon(:,:,:,1,1))),1),[])
%%
figure(1005); 
subplot(131); imshow(-270./log(squeeze(abs(P_recon(:,:,:,1,1)))./(squeeze(abs(P_recon(:,:,:,5,1)))+eps)),[1 500]); colormap jet
subplot(132); imshow(-270./log(squeeze(abs(P_recon(:,:,:,1,10)))./(squeeze(abs(P_recon(:,:,:,5,10)))+eps)),[1 500]); colormap jet
subplot(133); imshow(-270./log(squeeze(abs(P_recon(:,:,:,1,50)))./(squeeze(abs(P_recon(:,:,:,5,50)))+eps)),[1 500]); colormap jet
colorbar
%% 3D 
MR=MR2.Copy; K= sortArray(MR); size(K);
recon3d=bart('fft -i 7',double(K(:,:,:,1,1,1))); 
% recon3d=bart('rss -8',recon3d); 



%% TO DO 

 
 
x1=85
x2=115
x3=155
x4=180; 
y1=18
y2=20
y3=20
y4=20; 

figure(200);clf;subplot(221); hold on 
imshow(squeeze(abs(P_recon(:,:,1,6,1))),[]);
plot(y1,x1,'r+')
plot(y2,x2,'b+')
plot(y3,x3,'g+');
plot(y4,x4,'y+');
hold off
title('first echo, lowest TE');

subplot(223); hold on 
imshow(squeeze(abs(P_recon(:,:,1,1,end))),[]);
plot(y1,x1,'r+')
plot(y2,x2,'b+')
plot(y3,x3,'g+'); 
plot(y4,x4,'y+'); hold off
title('last echo, highest TE');
subplot(222)
hold on
plot([300,200,100,60,40,20,11],squeeze(abs(P_recon(x1,y1,1,:,1))),'r')
plot([300,200,100,60,40,20,11],squeeze(abs(P_recon(x2,y2,1,:,1))),'b')
plot([300,200,100,60,40,20,11],squeeze(abs(P_recon(x3,y3,1,:,1))),'g')
plot([300,200,100,60,40,20,11],squeeze(abs(P_recon(x4,y4,1,:,1))),'y')

hold off
title('T2prep dimension')
ylabel('intensity'); xlabel('TE')

subplot(224)
hold on
plot(squeeze(abs(P_recon(x1,y1,1,6,:))),'r')
plot(squeeze(abs(P_recon(x2,y2,1,6,:))),'b')
plot(squeeze(abs(P_recon(x3,y3,1,6,:))),'g')
plot(squeeze(abs(P_recon(x4,y4,1,6,:))),'y')

hold off
title('echo dimension')
ylabel('intensity'); xlabel('echo number')

%% T2 FITTING 
for i=1
T2PREP_data=gather((P_recon(:,:,1,:,i)));
size(T2PREP_data)
T2PREP=[300,200,100,70,50,30,10];
[T2_mono_all_allslice, T2_mono_all_rgb_allslice, rsquare_mono_all_allslice] = T2_fitting(T2PREP_data, T2PREP,0.02)
T2(:,:,i)=T2_mono_all_allslice; 
end
%%
figure(4000); imshow(median(T2,3),[10 400]); colormap parula 
colorbar; title('T2 map LRT')
figure(3001); hist(T2_mono_all_allslice((T2_mono_all_allslice>0)),200); title('T2 values')
 
 
 
 
 
 
 
 
 
 
 
 
