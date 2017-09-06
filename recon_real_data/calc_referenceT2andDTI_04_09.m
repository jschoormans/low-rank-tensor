clear; close all; clc
if ispc
    cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2017_08_15\')
else
    cd(['/home/',getenv('USER'),'/lood_storage/divi/Ima/parrec/Kerry/LRT_Data/2017_09_03']) %data is sorted in yesterday's folder
end


%% reference DTI fitting 
MR_DTI=MRecon('lr_04092017_1838208_21_2_wip_full-dtiV4.raw');
MR_DTI.Parameter.Parameter2Read.typ = 1;
disp('readdata')
MR_DTI.ReadData;
MR_DTI.RandomPhaseCorrection;
disp('corrections...')
% MR_DTI.RemoveOversampling;
MR_DTI.PDACorrection; %???
MR_DTI.DcOffsetCorrection;
MR_DTI.MeasPhaseCorrection;
disp('sortdata')
MR_DTI.SortData;
kspa = squeeze(MR_DTI.Data);
ima = bart('fft -i 7', kspa);


DTI_data = squeeze(bart('rss 8', ima));
DTI_data = circshift(DTI_data,21,2);

b = 500;
g = [0 0 0;...
    0.00000,  0.00000,  1.00000;...
    0.78225,  0.13995,  0.60704;...
    0.60142, -0.66611,  0.44112;...
    -0.29250,  0.74723,  0.59674;...
    -0.11759, -0.90381,  0.41147;...
    -0.92889,  0.28450,  0.23710;...
    0.51744,  0.80103,  0.30102;...
    -0.70473, -0.35950,  0.61166];


DTI_data = permute(abs(DTI_data),[2 3 1 4]);
[MD, FA, eigvec] = DTI_fitting(DTI_data, g, b);

mask = DTI_data(:,:,:,1)>1e5;
MD = MD.*mask;
FA = FA.*mask;
eigvec = bsxfun(@times,eigvec, mask );

figure(201); 
montage(permute(MD,[1 2 4 3]),'displayrange',[]); colormap jet; colorbar; title('MD');
figure(202); 
montage(permute(FA,[1 2 4 3]),'displayrange',[0 0.5]);  colorbar; title('FA'); colormap hot;
figure(203);montage(permute(squeeze(eigvec(:,:,:,:,1)),[1 2 4 3]));  colorbar; title('eigenvector #1');
figure(204);montage(permute(squeeze(eigvec(:,:,:,:,2)),[1 2 4 3]));  colorbar; title('eigenvector #2');
figure(205);montage(permute(squeeze(eigvec(:,:,:,:,3)),[1 2 4 3]));  colorbar; title('eigenvector #3');


%% reference T2 mapping
MR_t2prep = MRecon('lr_04092017_1831388_20_2_wip_full-t2prepV4.raw');
MR_t2prep.Parameter.Parameter2Read.typ = 1;
disp('readdata')
MR_t2prep.ReadData;
MR_t2prep.RandomPhaseCorrection;
disp('corrections...')
MR_t2prep.RemoveOversampling;
MR_t2prep.PDACorrection; %???
MR_t2prep.DcOffsetCorrection;
MR_t2prep.MeasPhaseCorrection;
disp('sortdata')
MR_t2prep.SortData;
kspa = squeeze(MR_t2prep.Data);
ima = bart('fft -i 7', kspa);


T2prep_data = squeeze(bart('rss 8', ima));
T2prep_data = circshift(T2prep_data,28,2);

T2prep_data = permute(flipdim(abs(T2prep_data), 4),[2 3 1 4]);

TEs = [10 20 40 60 80 100];
intensitycutoff = 1e5;
[T2_mono_all_allslice, T2_mono_all_rgb_allslice, rsquare_mono_all_allslice] = T2_fitting(T2prep_data, TEs, intensitycutoff);


