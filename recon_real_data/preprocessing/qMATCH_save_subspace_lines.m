clear all
addpath(genpath('/opt/amc/bart/')); vars;
cd('L:\basic\divi\Ima\parrec\Jasper\LRT\2017_09_19\qM_101945')
file='qm_19092017_2052467_15_2_wip_qmatch_invivo_allteV4.raw'
MR2=MRecon(file)
%%
MR2.Parameter.Parameter2Read.ky=[-2 -1 0 1 2].'
MR2.Parameter.Parameter2Read.kz=[-2 -1 0 1 2].'
MR2.Parameter.Parameter2Read.typ=1;
MR2.Parameter.Recon.ImmediateAveraging='No'
read_data_ix = (MR2.Parameter.Labels.Index.typ == 1);    

aver=repmat([0:199],[8 1]); 
aver=aver(:); 
aver=repmat(aver,[36 1]); 
aver=repmat(aver,[6 1]); 

aver_tem=MR2.Parameter.Labels.Index.aver;
aver_tem(read_data_ix) = aver(:);

MR2.Parameter.Labels.Index.aver=aver_tem; 
%%
MR22=MR2.Copy; 
MR22.Parameter.Parameter2Read.aver=[3:8:199].';
MR22.ReadData; 
MR22.RandomPhaseCorrection
MR22.RemoveOversampling
MR22.PDACorrection
MR22.DcOffsetCorrection
MR22.MeasPhaseCorrection
MR22.SortData
size(MR22.Data)

%%
subspacedata=squeeze(MR22.Data);
size(subspacedata)
subspacedata=subspacedata(:,106:110,17:21,:,:,:); 
imshow(abs(squeeze(subspacedata(100,:,:,1,1,1))),[0 1e-7]);
size(subspacedata)

save('subspacedata_qMATCH_09192017.mat','subspacedata')
%% do actual subspace
xrange=1:216
yrange=[3]
zrange=[1:5]
chans=[5]
TEs=1;
% subspacedataF=ifft(subspacedata,[],1);
% subspacedata_ch=bsxfun(@times,subspacedata,permute([-1,1,1,1,-1].',[2 1]));
% subspacedata_ch=bsxfun(@times,subspacedata_ch,permute([-1,1,1,1,-1].',[3 2 1]));
% size(subspacedata_ch)
size(subspacedata)

subspacedata2=reshape(subspacedata(xrange,yrange,zrange,chans,TEs,:),[length(xrange)*length(yrange)*length(zrange)*length(chans)*length(TEs),25]);
size(subspacedata2)

L4=4;
[nav_estimate_1,eigenvals_1]= subspace_estimator_multicoil(squeeze(subspacedata2),L4);

figure(1)
subplot(211)
plot(abs(nav_estimate_1))
subplot(212)
plot(eigenvals_1)