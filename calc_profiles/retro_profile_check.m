%% measured profile
cd('/home/qzhang/lood_storage/divi/Ima/parrec/Jasper/LRT/Low_Rank_2017_07_28');
MR=MRecon('lr_28072017_1150231_18_2_wip_vfa-dtiV4.raw')

read_idx = find(MR.Parameter.Labels.Index.typ == 1);
ky = MR.Parameter.Labels.Index.ky;
kz = MR.Parameter.Labels.Index.kz;
dyn = MR.Parameter.Labels.Index.dyn;

ky_read = ky(read_idx); ky_read = ky_read(1:15:end);
kz_read = kz(read_idx); kz_read = kz_read(1:15:end);
dyn_read = dyn(read_idx); dyn_read = dyn_read(1:15:end);

ky_range = (max(ky_read)- min(ky_read) + 1)
kz_range = (max(kz_read)- min(kz_read) + 1)
dyn_range = (max(dyn_read)- min(dyn_read) + 1)
kspa_frames = zeros(ky_range,kz_range,1,dyn_range);

for ii = 1:length(ky_read)
    ky_cor = ky_read(ii) - min(ky_read) + 1;
    kz_cor = kz_read(ii) - min(kz_read) + 1;
    dyn_cor = dyn_read(ii) - min(dyn_read) + 1;
    
    kspa_frames(ky_cor, kz_cor, 1, dyn_cor) = 1 + kspa_frames(ky_cor, kz_cor, 1, dyn_cor);
end

figure; montage(kspa_frames,'displayrange',[],'size',[6 9]); title('actual sampling pattern')

%% input profile
cd('/home/qzhang/lood_storage/divi/Ima/parrec/Kerry/CSENSE_profiles');
profile = dlmread('LRT_TSE_T2prep_64_64_9_6_r0_l1_bCtr5_sCtr2_us0.05.dat');

nr_profiles_per_dyn = 240;
ky_input_range = (max(profile(ii,1))- min(profile(ii,1)) + 1)
kz_input_range = (max(profile(ii,2))- min(profile(ii,2)) + 1)
dyn_input_range = ceil(length(profile)/nr_profiles_per_dyn)

profile_frames = zeros(ky_input_range,kz_input_range,1,dyn_input_range);

for ii = 1:length(profile)
    ky_input_cor = profile(ii,1) - min(profile(:,1)) + 1;
    kz_input_cor = profile(ii,2) - min(profile(:,2)) + 1;
    dyn_input_cor = ceil(ii / nr_profiles_per_dyn);
    
    profile_frames(ky_input_cor, kz_input_cor, 1, dyn_input_cor) = 1 + profile_frames(ky_input_cor, kz_input_cor, 1, dyn_input_cor);
end
figure; montage(profile_frames,'displayrange',[],'size',[6 9]); title('input sampling pattern')
    
    