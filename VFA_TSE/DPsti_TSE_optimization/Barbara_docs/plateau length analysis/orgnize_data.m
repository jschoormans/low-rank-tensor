figure;  set(0,'DefaultAxesFontSize', 16);

plot(squeeze(AllSignal(:,4,1)),squeeze(AllSignal(:,2,:)));
hold on

figure;
plot(squeeze(AllSignal(:,4,1)),squeeze(AllSignal(:,1,:)));

AllSignal=cat(3,plateau00, plateau01, plateau02, plateau03, plateau04, plateau05, plateau06, plateau07, plateau08, plateau09, ...
                plateau10, plateau11, plateau12, plateau13, plateau14, plateau15, plateau16, plateau17, plateau18, plateau19, ...
                plateau20, plateau21, plateau22, plateau23, plateau24, plateau25, plateau26);
            
            
%%
% k=0;            
eval(['plateau2' num2str(k),'=cat(2,[VarName2, VarName3, VarName4, VarName5, VarName6]);']); k=k+1
clear VarName1 VarName2 VarName3 VarName4 VarName5 VarName6 