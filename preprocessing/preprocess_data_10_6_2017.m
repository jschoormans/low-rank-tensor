% draft_code_to_load_and_preprocess_measured_data
cd('/home/jschoormans/lood_storage/divi/ima/parrec/jasper/LowRank_diffusion_apple_t2_2017_06_10')

file{1}='lo_10062017_1544080_2_2_wipdwitset2prep3V4.raw'
file{2}='lo_10062017_1618122_7_2_wipdwitseV4.raw'

%
addpath(genpath('/opt/amc/matlab/toolbox/MRecon-3.0.535')) % does not work
rmpath(genpath('/opt/amc/matlab/toolbox/MRecon'))
addpath(genpath('/opt/amc/matlab/toolbox/MRecon-3.0.519'))  % does not work 
rmpath(genpath('/opt/amc/matlab/toolbox/MRecon-3.0.506'))  % does not work
rmpath(genpath('/opt/amc/matlab/toolbox/MRecon-3.0.529')) % does not work AT ALL
rmpath(genpath('/opt/amc/matlab/toolbox/MRecon-3.0.482'))  % does not wrok
rmpath(genpath('/opt/amc/matlab/toolbox/MRecon-3.0.532'))  % does not work
rmpath(genpath('/opt/amc/matlab/toolbox/MRecon-3.0.515')) % does not work
rmpath(genpath('/opt/amc/matlab/toolbox/MRecon-3.0.523')) %other issue


clear MR
f=file{2}

MR=MRecon(f)
MR.Parameter.Parameter2Read.typ=1;
MR.Parameter.Recon.ACNrVirtualChannels=1
MR.Parameter.Recon.ArrayCompression='yes' %array compression so we dont have to deal with coil data for now
MR.Parameter.Recon.ImmediateAveraging='no'
MR.ReadData;
MR.DcOffsetCorrection;
MR.PDACorrection;
MR.RandomPhaseCorrection;
MR.MeasPhaseCorrection;
MR.SortData;
%%
MR.ShowData
MR.K2I;
MR.ShowData
%%

size(MR.Data)