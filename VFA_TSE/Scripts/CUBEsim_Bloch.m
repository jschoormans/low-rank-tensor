% Bloch Equation Simulation for CUBE

% This is build up using cells. Click in a cell and execute it with
% cntr + enter. The first cell initializes all values and here the
% simulation parameters can be set. The second cell performs the simulation
% and the other cells plot different results.

% The resulting vector M with size [3, total time(ms), number of spins]
% gives the magnetization for x, y and z at the end of each ms. 

% Thanks to Brian Hargreaves, Stanford,
% http://mrsrl.stanford.edu/~brian/bloch/
%
% Henk Smit, h.smit@erasmusmc.nl, 02-2013.
% -----------------------------------------

%%  ===== Initialize parameters ======
clear all

% ===== Set the parameters ======

% Tissue properties
nsp = 1000;         % number of spins with different random off-resonance freqs
T1 = 1000;        % ms.
T2 = 60;        % ms.
dfs= 1000;         % off-resonance scaling factor.

% Sequence Properties
TE = 4;          % ms, repetition time of FSPGR, only even numbers
ETL = 64;        % number of refocussing pulses
TRec = 30;        % time to wait after echotrains (so TR = TRec + ESP*ETL)
nET = 1;         % number of echo trains
invfa = 180;     % degree, flip angle inversion
spoilpart = 0;   % [0 1],  0=perfect spoiling, 1=no spoiling, spoiling after echo train before next inversion.
targsig = 0.25;   % target signal of CUBE, between 0 and 1, typically around 0.2

% ===== Get the Propagation Matrices ======

%calculate the off resonance frequencies [Hz]
df = dfs*(rand(1,nsp)-0.5);  

% df=-10:40/(nsp-1):30;

for fr=1:size(df,2) %free precession matrices for different dfs
    [K(fr,:,:),L(fr,:,:)] = freeprecess(1,T1,T2,df(fr));
end

excrot = yrot(90*2*pi/360);
% excrot = throt(90*2*pi/360,pi/2);
spoil = xyspoil(spoilpart);

% Get the flipangles for CUBE
angles=CUBEangles(T1,T2,TE,targsig,ETL);
anglesdeg=180*angles/pi;
TR = ETL*TE+TRec;  % time until next excitation
TEs= [];

%% ===== Simulate the sequence ======
% For large number of spins this can take a while

M = zeros(3,TR*nET,size(df,2));	% Keep track of magnetization at all time points.
Mrf = M; %to store the magn where the first ref RF acts on.
for ini=1:size(M,3)
	M(:,1,ini)=[0;0;1];	% Starting magnetization.
end

for sp=1:size(M,3)
    for i=0:nET-1
        if(i==0) %excitation pulse
            M(:,1,sp)=excrot*M(:,1,sp);
            M(:,1,sp) = squeeze(K(sp,:,:))*M(:,1,sp)+squeeze(L(sp,:,:))';
        else
            M(:,i*TR+1,sp) = spoil*M(:,i*TR,sp); %spoil by erasing xy magn before start next train
            M(:,i*TR+1,sp) = excrot*M(:,i*TR+1,sp);
            M(:,i*TR+1,sp) = squeeze(K(sp,:,:))*M(:,i*TR+1,sp)+squeeze(L(sp,:,:))';
        end

        for k=2:TR %k = time (ms)
            if(mod(mod(k+i*TR,TR),TE)<TE/2 && k<ETL*TE) %free precession before RF
                M(:,k+i*TR,sp) = squeeze(K(sp,:,:))*M(:,k+i*TR-1,sp)+squeeze(L(sp,:,:))';
            elseif(mod(mod(k+i*TR,TR),TE)==TE/2 && k<ETL*TE) %RF pulse
                TEs = [TEs k+i*TR+TE/2]; %store all acquisition times
                if mod((k+i*TR+TE/2),TR) ~= 0 %ugly fix for TRec=0
                    cuberot = xrot(angles(mod((k+i*TR+TE/2),TR)/TE)); %select the right angle(1:etl-1)
                else
                    cuberot = xrot(angles(ETL)); %select angle(etl)
                end
                M(:,k+i*TR,sp) = squeeze(K(sp,:,:))*M(:,k+i*TR-1,sp)+squeeze(L(sp,:,:))'; %free precession until RF
                
                if(k==TE/2 && i==0)
                    Mrf(:,k,sp)=M(:,k,sp);
                end
                
                M(:,k+i*TR,sp) = cuberot*M(:,k+i*TR,sp);  %RF
            elseif((mod(mod(k+i*TR,TR),TE)>TE/2 && k<=ETL*TE) || k==ETL*TE) %free precession after RF
                M(:,k+i*TR,sp) = squeeze(K(sp,:,:))*M(:,k+i*TR-1,sp)+squeeze(L(sp,:,:))';
            elseif(k>=ETL*TE) %free precession after echo train (recovery time)
                M(:,k+i*TR,sp) = squeeze(K(sp,:,:))*M(:,k+i*TR-1,sp)+squeeze(L(sp,:,:))';

            end
        end
    end
end



%% Plot net xy signal at acquisition times
clear i;

figure('position',[200 200 900 600])

%Calculate Complex Sum
Mcomb=squeeze(M(1,:,:)+i*M(2,:,:)); 

if size(M,3)>1
    meanMxyMag = abs(mean(Mcomb,2));
else
    meanMxyMag = abs(mean(Mcomb,1));
end

uTEs = unique(TEs);
plot(uTEs,meanMxyMag(uTEs),'r','LineWidth',3)
hold on
scatter(TEs(1:nET*ETL),zeros(1,nET*ETL),'.k')

legend(['Mean Mxy of ' num2str(nsp) ' spins'], 'Echo times');
xlabel('Time (ms)');
ylabel('Magnetization');
axis([0 max(uTEs) -0 1]);
grid on;

%% Plot mean Mz magn and Flip angles
% bug: requires #spins > 1

clear i;

figure('position',[200 200 900 600])

%Calculate Complex Sum
Mz=squeeze(M(3,:,:)); 
Mcomb=squeeze(M(1,:,:)+i*M(2,:,:)); 

if size(M,3)>1
    meanMzMag = abs(mean(Mz,2));
    meanMxyMag = abs(mean(Mcomb,2));
else
    meanMzMag = abs(mean(Mz,1));
    meanMxyMag = abs(mean(Mcomb,1));
end

uTEs = unique(TEs);
% plot(uTEs,meanMzMag(uTEs),'r','LineWidth',3)

hold on
% plot(meanMxyMag,'b','LineWidth',3)
plot(meanMzMag,'r','LineWidth',3)
% plot(uTEs,meanMxyMag(uTEs),'b','LineWidth',3)
scatter(TEs(1:nET*ETL),zeros(1,nET*ETL),'.k')
plot(TEs(1:nET*ETL),anglesdeg./360,'k','LineWidth',3)

% legend(['Mean Mxy of ' num2str(nsp) ' spins'],['Mean Mz of ' num2str(nsp) ' spins'], 'Echo times');
legend(['Mean Mz of ' num2str(nsp) ' spins'], 'Echo times','Flip angles/360');
xlabel('Time (ms)');
ylabel('Magnetization / Flip angle / 360 (degree)');
axis([-20 size(meanMzMag,1) -0 1]);
grid on;

%% Plot each spin in the xy plane right before the first rf pulse

% figure('position',[200 200 800 600])
% 
% for j=1:nsp
% %     g=plot([0 M(1,0+TE/2-1,j)],[0 M(2,0+TE/2-1,j)],'b','LineWidth',1.5);
%     g=plot([0 Mrf(1,0+TE/2,j)],[0 Mrf(2,0+TE/2,j)],'b','LineWidth',1.5);
%     hold on
%     cr=abs(randn(3,1));
%     set(g,'Color',cr-floor(cr))
% end
% 
% 
% 
% xlabel('Mx')
% ylabel('My')
% 
% axis([-1 1 -1 1])
% 
% title('All spins as vectors in the x-y plane')

%% Plot Mx, My, Mz for each off resonance in one figure
% For a large number of spins this takes a while. Plots each spin!

% clear i;
% 
% figure('position',[200 200 1000 700])
% time = [1:nET*TR];
% hold on
% 
% for pl=1:size(M,3)
% %     subplot(2,1,1)
%     plot(time,(M(1,:,pl)),'b-',time,(M(2,:,pl)),'r-',time,M(3,:,pl),'k-','LineWidth',2);
%     legend('M_x','M_y','M_z');
%     hold on
% end
    
