%% TO DO: 
% MOTION CORRECTION FOR LRT ALGORITHM 
% input k-space, output: corrected k-space 

% MAKE LOW RES RECONS (either row-wise, or variable time resolution)
% extract motion parameters (12 parameters would describe rigid motion
% completely?)
% correct k-space for this

%% draft code


% using brain image scanned at 24-1 

cd('L:\basic\divi\Ima\parrec\Jasper\LRT\Low_Rank_2018_01_24_VFA')

clear MR
MR=MRecon('lr_24012018_0731151_2_2_wip_brin-vfa-t2prep-lowrankV4.raw'); 
DTI=0;
% MR.Parameter.Parameter2Read.chan=[10; 11];
% MR.Parameter.Parameter2Read.ky=[-80:80].';

% MR.Parameter.Labels.Index.aver=(MR.Parameter.Labels.Index.rf);
% MR.Parameter.Parameter2Read.aver=[0:max(MR.Parameter.Labels.Index.aver)].'; 
MR.Parameter.Parameter2Read.typ = 1;
MR.Parameter.Recon.ImmediateAveraging='Yes';

disp('readdata')
MR.ReadData;
MR.RandomPhaseCorrection;

disp('corrections...')
MR.RemoveOversampling;
MR.PDACorrection; 
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;
MR.SortData


%% calculate sense maps 
K=permute(MR.Data,[1 2 3 4 6 5]);

sensmaps=bart('ecalib -m1',K);

imagine(sensmaps)
%% do low res recon 
bart('pics',K,sensmaps);

