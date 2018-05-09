function [trajGPU,w]=convert_radial_traj (traj)
% input bart traj [3,ky,kz] (2D)
[~,nspokes,res]=size(traj);

trajGPU=permute(traj,[3 2 1]); 
trajGPU=reshape(trajGPU,[nspokes*res,3]);
trajGPU=permute(trajGPU,[2 1]);
trajGPU=trajGPU./(res); % scale to -0.5 0.5 (DONT UNDERSTAND THIS YET) 


% for input in getradweights reshape in other way..
ku=squeeze(traj(1,:,:)+1i*traj(2,:,:)).';
w=getRadWeightsGA(ku);
w=col(w);

%output 
end