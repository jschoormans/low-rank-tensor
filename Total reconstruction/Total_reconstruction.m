%% CS-4Dflow-recon

% Editors (alphabetical order):
% l.m.gottwald@amc.uva.nl
% e.s.peper@amc.uva.nl
% v.q.pronk

% 17-October-2017

%% load raw data and perform compressed sensing or LRT reconstruction
clear all
close all
clc

% To check:
% Run 3D-sample + mrecon-steps afterwards and check if slice 7 is the same as just only a run on
% slice 7. Also check if location of combinecoils matters. Also try to do
% the way Eva and Lucas do it, uncombining coils after cs and recombining
% afterwards.
perform = 1;
retro_us = 0;
noCSorLRT = 0;
doLRT = 1; % Decides whether LRT and/or CS should be applied. 1 = yes
doCS = 1;
writeparrec = 0;
useperfectsensemaps = 0;
fixfoldover = 0;

% define folders and files
RP.name            = '17_04_2018_Flowphantomfully_y128z128c24_lrt_normalsensemaps20_nous_lambdamu0.8483opt_testtransformdomainrolePsi';
% RP.name            =econ_singleScript_hadam '06_02_2018_FlowPhantomFully_from2017_11_28_y64z64c24_3D_hadamardtransform_fullnavigator_withmasky1z1C20_rank6.6.4_2Dslice26';
% RP.data_dir        = '/home/barunderkamp/lood_storage/divi/Projects/4dflow_lrt/data/Ja_272623';
RP.data_dir        = '/home/barunderkamp/lood_storage/divi/Projects/4dflow_lrt/data/2017_11_28_FlowPhantomFully';
RP.data_target     = 'fl_28112017_1644108_5_2_wip_af_y64z64_c24_fullyV4';
% RP.data_target     = 'ja_08022018_1855474_5_2_wip_cs_r8_y128_z128_c17V4';
 
% % RP.data_dir        = '/home/barunderkamp/lood_storage/divi/Temp/Lukas/reconstruction/test-data/spiralshape';
% RP.data_target     = 'fl_10102017_1931244_11_2_wip_2dfully_30card_fov_80_nsa4_afprotocolV4';
% RP.data_senseref   = 'cs_30082017_1947257_1000_11_wip_senserefscanV4';
% RP.data_coilsurvey = 'cs_30082017_1946409_1000_8_wip_coilsurveyscanV4';


% initialize
init_toolbox;

% initialize MRecon object
mrecon = MRecon(strcat(RP.data_dir, filesep, RP.data_target,'.raw')); % load '.raw' or '.lab'
mrecon.Parameter.Recon.CoilCombination = 'pc';
% Let op: 07-02-2018 tijdelijk op Average gezet!!
% if retro_us;
    mrecon.Parameter.Cardiac.RetroHoleInterpolation='No';
% else 
%     RetroHoleInterpolation_Center
% end
%
mrecon.Parameter.Cardiac.Synchronization = 'Retrospective';
% mrecon.Parameter.Cardiac.RetroPhases = RP.Cardiac_RetroPhases; % define cardiac frames if needed
% mrecon.Parameter.Cardiac.RetroPhases = 30;
% define data to read
mrecon.Parameter.Parameter2Read.typ = 1;

%Added by Bobby 16-01-2018 as test, remove v_y and v_z encoding.
% mrecon.Parameter.Parameter2Read.extr1 = uint16([0;1]);
% mrecon.Parameter.Parameter2Read.extr1 = uint16([1]);
%
mrecon.Parameter.Parameter2Read.Update;

% load data
mrecon.ReadData;

% % apply retrospective undersampling
% mrecon.Data=mrecon.Data.*retro_mask;

% perform MRecon reconstruction steps
mrecon.RandomPhaseCorrection;
mrecon.RemoveOversampling;
mrecon.PDACorrection;
mrecon.DcOffsetCorrection;
mrecon.MeasPhaseCorrection;
mrecon.SortData;
mrecon.GridData; 


% Version from Jasper
% undersampling_mask=bart(' poisson -Y160 -Z160 -y5 -z5 -C15 -v -e');
% undersampling_mask = permute(undersampling_mask,[2 3 1]);
% mrecon.Data = bsxfun(@times,mrecon.Data,undersampling_mask);

% % Apply (retrospective) undersampling mask with different mask for each
% cardiac frame.

% if retro_us;
if ~perform;
% voeg toe: geen mask doen als fully sampled (perform, dus geen cs of lrt)
mrecon.Data = permute(mrecon.Data,[2,3,1,4:10]); % This to multiply mask by y-z planes.
seed = 1;
for cardframe = 1:size(mrecon.Data,6);
    for venc_dir = 1:size(mrecon.Data,10);
    mrecon.Data(:,:,:,:,:,cardframe,:,:,:,venc_dir) = bsxfun(@times,mrecon.Data(:,:,:,:,:,cardframe,:,:,:,venc_dir),permute(bart(sprintf(' poisson -Y64 -Z82 -y6 -z6 -C10 -v -e -s%d',seed)),[2 3 1]));
    seed = seed+1;
    end
end
mrecon.Data = permute(mrecon.Data,[3,1,2,4:10]); %... so remove if not taking x-slices.
end
% 
% save sampling mask
%

mask = mrecon.Data ~= 0;
% save_mat(dir_out(RP.data_dir), [RP.name '.mask'], 'mask', mask);

% perform MRecon reconstruction steps

mrecon.RingingFilter;
mrecon.ZeroFill;
mrecon.K2IM;
if noCSorLRT;
% mrecon.Data = mrecon.Data(7,:,:,:,:,:,:,:,:,:); % If taking only x-slice 7.
mrecon.K2IP;
mrecon.ConcomitantFieldCorrection;
mrecon.DivideFlowSegments; % after thorough search, found out this has to precede combinecoils.
mrecon.CombineCoils; % in CS done by BART
% tmp = bsxfun(@times,create_checkerboard([1,size(tmp,2),size(tmp,3)]),tmp); % apply checkerboard 





% mrecon.EPIPhaseCorrection ook erbij??
%TEMPORARY



% mrecon.CombineCoils; % in CS done by BART
mrecon.RemoveOversampling;

% Invert velocity direction (goes wrong somewhere)
%blablabla
%

save_mat( fullfile(dir_out(RP.data_dir)), [RP.name], '.Data', mrecon.Data);
        mrecon.WritePar([RP.name '.perform.par']);
        mrecon.WriteRec([RP.name '.perform.rec']);

else
tmp = mrecon.Data;
OversamplingFlag = mrecon.Parameter.ReconFlags.isoversampled;
RotationFlag = mrecon.Parameter.ReconFlags.isrotated;
end
if useperfectsensemaps;
perfect_sensemaps_path=fullfile('/home/barunderkamp/lood_storage/divi/Projects/4dflow_lrt/Bijmateriaal/Sensetests','fullysampledsensemaps_flowphantom.mat');
perfect_sensemaps_loaded = load(perfect_sensemaps_path);
perfect_sensemaps = perfect_sensemaps_loaded.perfect_sensemaps;
clear perfect_sensemaps_loaded
end
%% CS 
if doCS
%    tmpforCS.K2IM;
%    tmpforCS = tmpforCS.Data;
%     for i = [1:size(tmp,1)] % Loop over all 2D-slices
    for i = 26;

%     K_cs=tmpforCS(i,:,:,:,:,:,:,:,:,:,:,:);
    K_cs=tmp(i,:,:,:,:,:,:,:,:,:,:,:);

% De volgende stap zou correct een shift van -1,-1 moeten doen in image
% space en weer terug naar k-space moeten doen. Dit om de onverklaarbare
% shift te herstellen en fold-over te voorkomen.
% Door de fourier transforms krijgen de hoeken nu opeens wel waarden niet-nul \n')
% misschien beter om dat eerst nul te maken. zelfde kwestie ongeveer als kerry zei \n')
% maar daar kijken we later wel naar. bekijk eerst resultaat maar even \n')

%part below was attempted but failed.

% fprintf('check hier even of de mask er nog wel is alvorens de shift')
% fprintf('oke ik weet al waarom dit niet werkt denk ik. omdat de mask niet \n')
% fprintf(' meer een mask is omdat die wordt verspreid door de fourier transform. hoe dit op te lossen?')
% K_cs=fft(((fft((circshift((sqrt(size(K_cs,3)).*(ifft((sqrt(size(K_cs,2)).*(ifft(K_cs,[],2))),[],3))),[0,-1,-1])),[],2))./sqrt(size(K_cs,2))),[],3)./sqrt(size(K_cs,3));
% 
% fprintf('test of checkerboard nog goed werkt na shift want er gebeurt iets fouts \n')
% fprintf('ben ook checkende of het ligt aan lambda');
% apply checkerboard pattern 
%attempts to shift image 1 pixel diagonally.
if fixfoldover;
fprintf('beter om checkerboard-style te vermenigvuldigen qua phase')
K_cs = bsxfun(@times,create_checkerboard_2(),K_cs);
end
%
K_cs = bsxfun(@times,create_checkerboard([1,size(K_cs,2),size(K_cs,3)]),K_cs); % 13-12-2017, changed mrecon.Data into K_cs in 'size'  


% apply coil compresseion to 8 coils
% usage: cc [-p d] [-M] [-r ...] [-A] [-S ...] [-G ...] [-E ...] <kspace> <coeff>|<proj_kspace>
% cmdcc = 'cc -p 8 -E';        
% [T,tmpforCS] = evalc('bart(cmdcc,tmpforCS)');

if ~useperfectsensemaps;
% estimate sensitivity maps
% usage: ecalib [-t f] [-c f] [-k ...] [-r ...] [-m d] [-S] [-W] [-I] [-1] [-v f] [-a] or: caldir [cal_size]  
cmdsens = 'ecalib -r20 -m1 -S -I';
% [T,sensemapforCS] = evalc('bart(cmdsens, sum(K_cs(:,:,:,:,1,:,1,1,1,1),6)./sum(K_cs(:,:,:,:,1,:,1,1,1,1)~=0,6) )');
[T,sensemapforCS] = evalc('bart(cmdsens, sum(sum(K_cs(:,:,:,:,1,:,1,1,1,:),6),10)./sum(sum(K_cs(:,:,:,:,1,:,1,1,1,:)~=0,6),10) )');
else
    sensemapforCS = perfect_sensemaps(:,:,:,:,i);
end

    % do compressed sensing reconstruction
    % usage: pics [-l ...] [-r f] [-c] [-s f] [-i d] [-t <string>] [-n] [-g] [-p <string>] [-I ...] [-b d] [-e] [-W <string>] [-d d] [-u f] [-C d] [-f f] [-m ...] [-w f] [-S] [-B d] [-K] <kspace> <sensitivities> <output>
%     cmdpics = 'pics -R W:7:0:0.0251 -R T:1024:0:0.1585 -i 50 -S -d 5'; 
    cmdpics = 'pics -R T:1024:0:0.0316 -i 50 -S -d 5'; 

    K_cs=permute(K_cs, [1:5 7:10 11 6]); % permute 
    [T,K_cs] = evalc('bart(cmdpics,K_cs,sensemapforCS)'); 
    
    % finalize Data_onlycs
    Data_onlycs = permute(single(K_cs), [1:5 11 6:9 10]); % permute back  
    
    % stack all slices to 3d dataset
Data_onlycs = squeeze(Data_onlycs);
p = ndims(Data_onlycs)+1;

if i == 26;
% if i == 1
    Data_onlycsstack = Data_onlycs;
else
    Data_onlycsstack = cat(p,Data_onlycsstack,Data_onlycs);
end

    end
    
% Permute x-stacks to traditional 1st index
Data_onlycsstack = permute(Data_onlycsstack, [p,1:p-1]);

    % save Data_onlycsstack after complex division (fix this) -> fixed?
    mrecon.Data = permute(Data_onlycsstack,[1:3,6,7,4,8,9,10,5]);
    mrecon.Parameter.ReconFlags.isimspace = [1,1,1];
    mrecon.Parameter.ReconFlags.isoversampled = OversamplingFlag;
    mrecon.Parameter.ReconFlags.isconcomcorrected = 0;
    mrecon.ConcomitantFieldCorrection;
    mrecon.DivideFlowSegments;
    mrecon.RemoveOversampling;   
   
%     Data_cs_compldiv=permute(mrecon.Data,[1:3,6,10,4,5,7,8,9]); % Squeeze but keep x-dimension (size 1)
Data_cs_compldiv = mrecon.Data;
%     for j = [1:3];% alternative method
%         Data_cs_compldiv(:,:,:,:,j) = Data_onlycsstack(:,:,:,:,1)./Data_onlycsstack(:,:,:,:,j+1);
%     end

    Data_cs_compldiv_velocityinverted = ((Data_cs_compldiv).^(-1)).*(abs(Data_cs_compldiv).^2);
    
    mrecon.Data = Data_cs_compldiv_velocityinverted;
    %temporary save 09-04-2018
    save_mat( '/home/barunderkamp/lood_storage/divi/Projects/4dflow_lrt/data/2017_11_28_FlowPhantomFully/recon_out/CS tests amount of iters and wavelike artefact', [RP.name '.DataforCS_phasedif'], 'Data', Data_cs_compldiv_velocityinverted);
    %
    save_mat( fullfile(dir_out(RP.data_dir)), [RP.name '.DataforCS_nophasedif'], 'Data', Data_onlycsstack); % this one still has oversampling!
    save_mat( fullfile(dir_out(RP.data_dir)), [RP.name '.DataforCS_phasedif'], 'Data', Data_cs_compldiv_velocityinverted);
    if writeparrec;
    % Save parrec
    mrecon.WritePar([RP.name '.CS.par']);
    mrecon.WriteRec([RP.name '.CS.rec']);
    end
end
%
%% LRT
if doLRT
%     tmpforLRT.K2IM; % If 2D-dataset, remove this line, for-loop and concatenation.
%     tmpforLRT=tmpforLRT.Data;
%     for i = [1:size(tmp,1)] % Loop over all 2D-slices;
    for i = [26];
        close all
%     K=tmpforLRT(i,:,:,:,:,:,:,:,:,:,:,:);
    K=tmp(i,:,:,:,:,:,:,:,:,:,:,:);
    
    % oke dit werkt niet voor LRT, want dan wordt de mask verspreid tot
    % alle punten.
    fprintf('shift werkt niet voor LRT want mask wordt verspreid tot alle punten')
% K_cs=fft(((fft((circshift((sqrt(size(K_cs,3)).*(ifft((sqrt(size(K_cs,2)).*(ifft(K_cs,[],2))),[],3))),[0,-1,-1])),[],2))./sqrt(size(K_cs,2))),[],3)./sqrt(size(K_cs,3));

    
% remove stupid checkerboard pattern
% che=create_checkerboard([160,160,1]);
che=create_checkerboard([1,size(K,2),size(K,3)]); % 

K=bsxfun(@times,K,che);

% % apply coil compresseion to 8 coils
% % usage: cc [-p d] [-M] [-r ...] [-A] [-S ...] [-G ...] [-E ...] <kspace> <coeff>|<proj_kspace>
% cmdcc = 'cc -p 8 -E';        
% [T,K] = evalc('bart(cmdcc,K)');

% EVA'S WAY (AS IN CS):
% % estimate sensitivity maps
% % usage: ecalib [-t f] [-c f] [-k ...] [-r ...] [-m d] [-S] [-W] [-I] [-1] [-v f] [-a] or: caldir [cal_size]  
% cmdsens = 'ecalib -r20 -m1 -S -I';
% [T,sensemapforLRT] = evalc('bart(cmdsens, sum(K(:,:,:,:,1,:,1,1,1,1),6)./sum(K(:,:,:,:,1,:,1,1,1,1)~=0,6) )');

%ORIGINAL WAY (AS IN FLOW_RECON_2):
kspace = squeeze(K);
if ~useperfectsensemaps;
a = sum(sum(kspace,4),5)./sum(sum(kspace~=0,4),5);
sens=bart('ecalib -r 20 -m1 -S -I',permute(a,[4 1 2 3]));
fprintf('neem wel dezelfde -r calibration region grootte als bij CS!')
else sens=perfect_sensemaps(:,:,:,:,i);
end
[P_recon,goldstandardscaled,C,G1,G2,G3] = LRT_recon_step1(kspace,sens);

% stack all slices to 3d dataset
p = ndims(P_recon)+1;

if i == 26
% if i == 1
    P_reconstack = P_recon;
else
    P_reconstack = cat(p,P_reconstack,P_recon);
end
%
    end
    % Permute x-stacks to traditional 1st index

    P_reconstack = permute(P_reconstack, [p,1:p-1]);
    %

    % save P_reconstack (after doing complex division) Fix this. -> fixed?
    
    mrecon.Data = permute(P_reconstack,[1:3,6,7,4,8,9,10,5]);
    mrecon.Parameter.ReconFlags.isimspace = [1,1,1];
    mrecon.Parameter.ReconFlags.isoversampled = OversamplingFlag;
    mrecon.Parameter.ReconFlags.isconcomcorrected = 0;
    mrecon.ConcomitantFieldCorrection;
    mrecon.DivideFlowSegments;
    mrecon.RemoveOversampling;

%     P_reconcompldiv=permute(mrecon.Data,[1:3,6,10,4,5,7,8,9]);
P_reconcompldiv = mrecon.Data;
%     for j = [1:3]; % Alternative method
%         P_reconcompldiv(:,:,:,:,j) = P_reconstack(:,:,:,:,1)./P_reconstack(:,:,:,:,j+1);
%     end

    P_reconcompldiv_velocityinverted = ((P_reconcompldiv).^(-1)).*(abs(P_reconcompldiv).^2);
    
    mrecon.Data = P_reconcompldiv_velocityinverted;


%     save_mat( fullfile(dir_out(RP.data_dir)), [RP.name '.DataforLRT_nophasedif'], 'Data', P_reconstack); % RECON-BART RESULT FROM FLOW_RECON_2 NOT SAVED YET! | This one still has oversampling!
    save_mat( fullfile(dir_out(RP.data_dir)), [RP.name '.DataforLRT_phasedif'], 'Data', P_reconcompldiv_velocityinverted);

    if writeparrec == 1;
        % Save parrec
        mrecon.WritePar([RP.name '.LRT.par']);
        mrecon.WriteRec([RP.name '.LRT.rec']);
    end

end
%% Save parrec
% mrecon.Data=data;
%         mrecon.Data = permute(mrecon.Data,[1 2 3 6 7 4 8 9 10 5]); % permute for par rec export
        mrecon.WritePar([RP.name '.perform.par']);
        mrecon.WriteRec([RP.name '.perform.rec']);
        %%
%% continue here if you want to load the measured sensitiviy maps and use them for the final mrecon steps
% NOT FIXED YET TO HANDLE LRT AND/OR CS IN THIS COMBINED SCRIPT.

% load measured sensitiviy maps
S = MRsense(strcat(RP.data_dir, RP.data_senseref,'.raw'), strcat(RP.data_dir, RP.data_target,'.raw'), strcat(RP.data_dir, RP.data_coilsurvey,'.raw'));
S.OutputSizeReformated = [mrecon.Parameter.Encoding.XReconRes, ...
                          mrecon.Parameter.Encoding.YReconRes,...
                          mrecon.Parameter.Encoding.ZReconRes];
S.OutputSizeSensitivity = S.OutputSizeReformated;
S.Mask = 1;
S.Smooth = 1;
S.Extrapolate = 1;
S.Perform;
mrecon.Parameter.Recon.Sensitivities = S;
mrecon.Parameter.Recon.SENSERegStrength = 0;
sensemap = mrecon.Parameter.Recon.Sensitivities.Sensitivity;
clear S

% performe fft and uncombine coils
tmp = bart('fft -u 7',bart('fmac -s 16',tmp,sensemap)); % first fmac then fft
tmp = permute(single(tmp), [1:5 11 6:9 10]); % permute back
tmp = bsxfun(@times,create_checkerboard([1,size(tmp,2),size(tmp,3)]),tmp); % apply checkerboard 
mrecon.Data = tmp;

% perform MRecon reconstruction steps
mrecon.K2IM;
mrecon.EPIPhaseCorrection; 
mrecon.K2IP;
mrecon.GridderNormalization; 
mrecon.SENSEUnfold;
mrecon.PartialFourier;
mrecon.ConcomitantFieldCorrection;
mrecon.DivideFlowSegments;
mrecon.CombineCoils;
mrecon.Average;
mrecon.GeometryCorrection;
mrecon.RemoveOversampling;
mrecon.ZeroFill;
mrecon.FlowPhaseCorrection;
mrecon.RotateImage;
Data_allsteps = mrecon.Data; 

% save Data_allsteps
save_mat( fullfile(dir_out(RP.data_dir)), [RP.name '.Data'], 'Data', Data_allsteps  );



