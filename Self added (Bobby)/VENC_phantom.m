function Pvenc = VENC_phantom()
close all
N=128;
for nphantom=1:4
    E(1)=0.4+0.6*rand; %intensity
    E(2)= 0.1+rand*0.2; %length
    E(3)= 0.1+rand*0.2; %width
    E(4)=-0.5+(rand); %x-coord of middle 
    E(5)=-0.5+(rand); %y-coord of middle
    E(6)= rand*360; %angle   
    Pvenc{nphantom} = phantom(E,N);
    
    
%     TSE_Mxy_modulation(nphantom,:) = CUBE_Mxy_calculation_lite(T1vals(nphantom)*1e3, T2vals(nphantom)*1e3, TEes*1e3, ETL); % this function take time unit of ms instead of s
end