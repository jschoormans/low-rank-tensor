addpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations\tensor\low-rank-tensor'))

params=initprofile_params;
params.undersampling=240;
mask= make_profile_function(params);


%% order and save

makefilename = @(x) [x,num2str(params.ky),'_',num2str(params.kz),'_',num2str(params.nDim1),'_',num2str(params.nDim2),...
        '_r',num2str(params.radialflag),'_l',num2str(params.linearflag),'_bCtr',num2str(params.bigctrsize),...
        '(',num2str(params.dim1_bigctr),',',num2str(params.dim2_bigctr),')_sCtr',num2str(params.smallctrsize),'_us',num2str(params.undersampling)];
    

if ~params.DTIflag
    profile_order=profile_ordering(mask,params.radialflag,params.linearflag,params.visualize);
    filename=makefilename('LRT_VFA_T2p_'); 
else
    shot_per_frame = params.nr_points./params.ETL;
    assert((round(shot_per_frame)==shot_per_frame),'Profiles in every dynamic must be acquired in integer no. of shots! Change nr_points or ETL!');
    profile_order=profile_ordering_DTI_t2prep(mask,radialflag,linearflag,visualize,shot_per_frame);
    filename=makefilename('LRT_DTI_T2p_'); 
end

 
if ispc()
     cd('L:\basic\divi\Ima\parrec\Jasper\profiles_LRT')
else
     cd(['/home/',getenv('USER'),'/lood_storage/divi/Ima/parrec/Jasper/profiles_LRT']);
end
savemask_LRT(profile_order,filename,visualize)
fprintf('Saved as %s \n',filename)
fprintf('in folder: %s \n',pwd)
fprintf('Finished! \n \n \n ')
