% make synthetic radial (golden-angle) data 
% goal: 5D radial DCE\

addpath(genpath('C:\Users\jschoormans\gpuNUFFT-master'));
addpath(genpath('L:\basic\divi\Projects\cosart\tensor\low-rank-tensor'))
vars
%%
res=128
nspokes=256
nch=8
traj=bart(['traj -x',num2str(res),' -y',num2str(nspokes),' -r -G']);
ksp=bart(['phantom -k -s',num2str(nch),' -t'],traj);
sens=bart(['phantom -S',num2str(nch),' -x',num2str(res)]);

%% try to fit bart coordinates into GPUNUFFT coords 

osf =1;
wg = 3; sw = 8;

[trajGPU,w]=convert_radial_traj(traj);

FT = gpuNUFFT(trajGPU,w,osf,wg,sw,[res,res,1],sens,true);

im=FT'*reshape(ksp,[res*nspokes,nch]);
figure(1); imshow(abs(im),[])
%% make phantom images
res=128
nframes=50;
period=10;

I=DCE_cardiac_phantom(res,nframes,period);

% sort in radial trajectories (not you can also sort in 'true time',
% instead of 'binnned timebins'--> more exact
I_sorted=reshape(I,[res res 1 1 period nframes/period]);
%% kspace trajectory
spokesperframe=11; 
traj=bart(['traj -x',num2str(res),' -y',num2str(nframes*spokesperframe),' -r -G']);

traj_sorted=permute(traj,[1 2 3]); 
traj_sorted=reshape(traj_sorted,[3 res spokesperframe period nframes/period]); % not sure about ordering
traj_sorted=permute(traj_sorted,[1 3 2 4 5]); %[3 nspokes res, period, timeframes]

%% bin data
size(I_sorted); %    64    64     1     1    10    25
size(traj_sorted); %      3    11    64    10    25


clear trajGPU w
for ii=1: period
    for jj=1: nframes/period
        [trajGPU(:,:,ii,jj),w(:,ii,jj)]=convert_radial_traj(traj_sorted(:,:,:,ii,jj));
    end
end

osf=2
MDFT = MDgpuNUFFT(trajGPU,w,osf,wg,sw,[res,res,1],[],true);
radialksp=MDFT*I_sorted ;
temp2=MDFT'*(radialksp.* sqrt(permute(w,[1 4 2 3])));
imshow(abs(temp2(:,:,1,1,1,1)),[])


%%
save('synthetic_radial.mat','trajGPU','radialksp','w','I_sorted','res','nframes','period')
