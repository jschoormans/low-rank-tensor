%%
clear; clc; close all
vars_kerry;
%%
clear; 
current_mat_file = 'sc_27.mat'
raw_data_fn = ['/home/qzhang/lood_storage/divi/Ima/parrec/Jasper/Low_Rank_2017_06_12/',...
    'lo_12062017_2051218_27_2_wipdwitset230V4.raw'];

MR_data_raw = MRecon(raw_data_fn);

MR_data = MR_data_raw.Copy;
MR_data.Parameter.Parameter2Read.typ = 1;
MR_data.Parameter.Parameter2Read.mix = 0;  
MR_data.Parameter.Recon.ImmediateAveraging = 'Yes';
MR_data.ReadData;
MR_data.RandomPhaseCorrection;
% MR_data.RemoveOversampling;
MR_data.PDACorrection;
% MR_data.DcOffsetCorrection;
MR_data.MeasPhaseCorrection;



%==============Continue with Image Kspa data extracting
MR_data.SortData;
% >>>>>>>exported kspace data from here
k_spa_data = squeeze(double(MR_data.Data));

MR_data.GridData;
MR_data.RingingFilter;
MR_data.ZeroFill;

MR_data.K2IM;
MR_data.EPIPhaseCorrection;

MR_data.K2IP;
MR_data.GridderNormalization;

%--------------Calculate SENSE object-------------
% sense_ref = 'dp_24052017_1650436_1000_11_wipsenserefscanexperiment1clearV4.raw';
% MR2 = raw_data_fn;
% coil_survey = 'dp_24052017_1649131_1000_8_wipcoilsurveyscanexperiment1V4.raw';
% MR_sense = MRsense(sense_ref, MR2, coil_survey);
% MR_sense.Mask = 1;
% MR_sense.MatchTargetSize = 1;
% MR_sense.Perform;
% MR_DPstiTSE.Parameter.Recon.Sensitivities = MR_sense;
%----------------------end-----------------------


% MR_data.SENSEUnfold;  %no sense here

MR_data.ConcomitantFieldCorrection;
MR_data.DivideFlowSegments;
MR_data.Parameter.Recon.CoilCombination = 'pc';
MR_data.CombineCoils;
MR_data.ZeroFill;
MR_data.FlowPhaseCorrection;
MR_data.RotateImage;

MR_data.ShowData;
ima_data_cc = squeeze(double(MR_data.Data));
size(ima_data_cc)
% eval(['save ', current_mat_file, 'ima_data -v7.3']);

