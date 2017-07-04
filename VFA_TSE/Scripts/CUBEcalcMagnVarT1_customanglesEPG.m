%Script to show signal evolution of different tissues during CUBE. CUBE flip angles are
%calculated for one T1 and T2, script visualizes how other tissues respond.
%
%Based on the paper by Busse 2006
%Called functions:
%function [Mss] = CUBEmagn(T1,T2,TE,Angles,etl,StMagn)
%function [Angles,Mss] = CUBEangles(T1,T2,TE,Star,etl)

clear all
close all

%parameters to calculate CUBE flip angles
T1 = 1000;
T2 = 80;
dTE = 3;
Star = 0.20;
DegAngles = load('D:\Bram\MATLAB\CUBE Simulations\VFA_output_scanner\flips_PDVar.csv');
DegAngles = DegAngles(2:end)';
%etl = 61;
etl = length(DegAngles);
B1factor = 1.0;
TEprep = 1;

%parameters of simulated tissues 
%T1s and T2s need to have the same length!
T1s = [200 400 600 800 1000];
T2s = [50 50 50 50 50];

%Calculate Mxy for 5 different tissues
figure(1)
for j=1:size(T1s,2)
    % Get flip angles
    Angles = DegAngles.*(pi/180);
    % Plot flip angles 
    subplot(1,2,1)
    axis([0 etl 0 190])
    % legend('Flip angles')
    xlabel('Echo')
    ylabel('Flip Angle (degree)')
    % title('Flip angle evolution for different T2s')
    plot(1:size(DegAngles,2),DegAngles,'Color',j*[0.03 0.15 0.1],'LineWidth',4)
    hold on
    
    % Get magnetization
    [FM,ZM] = CUBEmagn(T1s(j),T2s(j),dTE,Angles,etl,exp(-TEprep/T2s(j)));
    sig(j,:)=FM(round(size(FM,1)/2),:);
   
    % Plot magnetization
    subplot(1,2,2)
    xlabel('Echo')
    ylabel('Magnetization')
    axis([0 75 0 0.8])
    plot(1:size(Angles,2),sig(j,:),'Color',j*[0.03 0.15 0.1],'LineWidth',3.5)
    hold on
    %scatter(1:size(Angles,2),sig(j,:),'.b','SizeData',350) 
end

subplot(1,2,2)
legend(strcat('T1=',num2str(T1s(1)),'  T2=',num2str(T2s(1))),...
        strcat('T1=',num2str(T1s(2)),'  T2=',num2str(T2s(2))),...
        strcat('T1=',num2str(T1s(3)),'  T2=',num2str(T2s(3))),...
        strcat('T1=',num2str(T1s(4)),'  T2=',num2str(T2s(4))),...
        strcat('T1=',num2str(T1s(5)),'  T2=',num2str(T2s(5))));

