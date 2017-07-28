%calculate scan profile for LRT scan 
%where dimension 1 is the TSE echo number and dimension 2 is T2prep

%%%%%%%%%%%% PARAMETERS TO CHANGE%%%%%%%%%%%%%%%%%%%%%%%
nDim1=20; % TSE dimensions
nDim2=5; % T2-prep 
ky=64; 
kz=64; 

bigctrsize=5;
smallctrsize=2;
undersampling=0.05; %excluding centers
% nr_points =240;

% waiting_time=(521-127)e-3; 
% TR=5.21e-3         %TR in ms; 
% TR_shot=nDim1*TR+waiting_time; 
TR_shot=521e-3

MC_maxiter=10000; 

visualize=1
radialflag=0 %radial/linear
linearflag=1; % 0 vertical ordering/ 1 horizontal ordering;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculating some params...
nr_centerpoints=(2*bigctrsize+1)^2; %number of k-points in the center squares; 
nr_points=ceil(undersampling*ky*kz) 
assert(nr_points>=nr_centerpoints,'fully sampled centers too big relative to undersampling')
nshots=nDim2*nr_points; 
total_time= TR_shot*nshots; %total time in seconds; 
fprintf('number of shots: %d, TSE number: %d, total time: %d seconds \n',nshots,nDim1,round(total_time))

%% add random point to every independent k-space
% performs a Monte Carlo simulation s.t. all have equal # of points
mask=zeros(ky,kz,nDim1,nDim2); 
fprintf('Starting Monte Carlo simulation \n')
for dim1=1:nDim1;
    for dim2=1:nDim2;
        fprintf('Dim 1: %d, Dim 2: %d |',dim1,dim2)
        m=zeros(ky,kz); MC_niter=0;
        if dim1==1 || dim2==1;
            ctrsize=bigctrsize;
            nr_centerpoints=(2*ctrsize+1)^2; %number of k-points in the center squares; 
        else
            ctrsize=smallctrsize;
            nr_centerpoints=(2*ctrsize+1)^2; %number of k-points in the center squares; 

        end
            while sum(m(:))~=nr_points
                m=rand(size(m))>(1-((nr_points-nr_centerpoints)/(ky*kz)));
                m=addCtr(m,ctrsize);

                MC_niter=MC_niter+1;
                if MC_niter>MC_maxiter
                    error('too many MC iterations - check settings')
                end
                
            end
        fprintf(' Monte Carlo iterations: %d \n',MC_niter)
        mask(:,:,dim1,dim2)=m;
    end
end
clear m
if visualize;   figure(1); clf; imshow(reshape(permute(mask(:,:,1:2,1:2),[1 3 2 4]),[ky*2,kz*2])); end

%% order and save
% profile_order=profile_ordering(mask,radialflag,linearflag,visualize);
profile_order=profile_ordering_DTI_t2prep(mask,radialflag,linearflag,visualize);

filename=['LRT_TSE_T2prep_',num2str(ky),'_',num2str(kz),'_',num2str(nDim1),'_',num2str(nDim2),...
     '_r',num2str(radialflag),'_l',num2str(linearflag),'_bCtr',num2str(bigctrsize),'_sCtr',num2str(smallctrsize),'_us',num2str(undersampling)]
 
if ispc()
     cd('L:\basic\divi\Ima\parrec\Jasper\profiles_LRT')
 else
end
savemask_LRT(profile_order,filename,visualize)

