%Script to show flip angle evolution with CUBE for different T2s
%Based on Extended Phase Graph algorithm as described in paper by Busse 2006

%Called functions:
%function [Mss] = CUBEmagn(T1,T2,TE,Angles,etl,StMagn)
%function [Angles,Mss] = CUBEangles(T1,T2,TE,Star,etl)

clear all

%parameters to calculate CUBE flip angles
T1 = 1000;
T2 = 50;
Star = 0.3;
TE = 5;
etl = 30;

a=figure('Position',[100 100 600 400]);

for j=1:size(T2,2)
    [Angles]=CUBEangles(T1,T2(j),TE,Star,etl);
    DegAngles=Angles.*180/(pi);
    plot(1:size(DegAngles,2),DegAngles,'Color',j*[0.03 0.15 0.1],'LineWidth',4)
    hold on
    
end

axis([0 etl 0 190])
% legend('Flip angles')
xlabel('Echo')
ylabel('Flip Angle (degree)')
% title('Flip angle evolution for different T2s')

set(gca,'fontsize',9)
name = ['Angles.tif'];
box off
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 2.9])



