% Bloch Equation Simulation for CUBE

% Results in a plot of the Mx, My, Mz, mean Mxy, and phase

% All relevant variables can be set in the first lines. 

% The first variable, nsp, sets the number of different off resonances. 
% for a clear subplot 1 and 3 make nsp small
% for a clear subplot 2 make nsp large
%
% Thanks to Brian Hargreaves, Stanford,
% http://mrsrl.stanford.edu/~brian/bloch/
%
% Henk Smit, h.smit@erasmusmc.nl, 09-2012.
% -----------------------------------------

clear all; close all
set(0,'DefaultAxesFontSize', 16);

 TErep=1;

    % ===== Set the parameters ======

    % Tissue properties
    nsp = 20;         % number of spins with different random off-resonance freqs
    
    %carotid vw
    T1 = 1000;        % ms.
    T2 = 50;        % ms.
%     
%     %nerve 
%     T1 = 1080;        % ms.
%     T2 = 78;        % ms.
    
    %muscle neck
%     T1 = 900;        % ms.
%     T2 = 35;        % ms.
%     
%     %fat
%     T1 = 382;        % ms.
%     T2 = 133;        % ms.

%     %white matter
%     T1 = 1010;
%     T2 = 90; 
        
    
    df0 = 0;         % df0=1 will only display the results for one on resonance spin


    % Sequence Properties
    TE = 7.3;          % ms, repetition time of FSPGR
    ETL =30;        % number of inversions
    TRec = 10;      % time to wait after echotrains
    nET = 1;         % number of echo trains
    invfa = 180;     % degree, flip angle inversion
    spoilpart = 0;   % [0 1],  0=perfect spoiling, 1=no spoiling, spoiling after echo train before next inversion.

    % T2 prep properties // 90x - 180y 180y (...) - -90x
    T2ETL = 4;      % number of refocussing pulses
    T2int = 6;      % time between refocussing pulses
    b1 =  180/180;       % factor multiplied with refocussing pulses (1 = 180 degree)


    %  ===== Get the Propagation Matrices ======


    %calculate the off resonance frequencies [Hz]
    df = 1300*(rand(1,nsp)-0.5);  

    if df0
        df=0;
    end

    if size(df,2)==1
        df=[df df];
    end

    for fr=1:size(df,2) %free precession matrices for different dfs
        [K(fr,:,:),L(fr,:,:)] = freeprecess(0.01,T1,T2,df(fr));
    end

    invrot = xrot(invfa*2*pi/360);
    excrot = yrot(90*2*pi/360);


    spoil = xyspoil(spoilpart);

    % Get the flipangles for CUBE
%     angles=CUBEangles(T1,T2,TE,0.24,ETL);
%     angles=CUBEangles(1200,80,TE,0.2,ETL);
%     angles = b1 *angles;
% 
%     angdeg=angles*180/pi;
%     plot(1:length(angdeg),angdeg)

    %input angle
%    angg=[180 156 127 120 120 120 120 121.6 123.2 124.8 126.4 128 129.6 131.2 132.8 134.4 136 137.6 139.2 140.8 142.4 144 145.6 147.2 148.8 150.4 152 153.6 155.2 156.8 168.4 180] %TSE 32 Low-High
    %angg=[180 137.25 84.5 64.31 58.13 58.13 63.1 71.4 83.42 101.52 120 121.9 123.81 125.71 127.62 129.52 131.43 133.33 135.24 137.14 139.05 140.95 142.86 144.76 146.67 148.57 150.48 152.38 154.29 156.19 168.1 180] %TSE 30 Linear
    %angg = [180 154.13 117.63 106.25 106.25 106.25 120 122.86 125.71 128.57 131.43 134.29 137.14 140 142.86 145.71 148.57 151.43 154.29 167.14 180];  %TSE 20
%     angg = [180 139.93 74 49.89 43.67 35*ones(1,112), 40.67 49.89 74 139.39 180 ];
%     angg = [137.25 84.5 64.31 58.13 58.13 63.66 73 87.18 111.01 120 121.21 122.42 123.64 124.85 126.06 127.27 128.48 129.7 130.91 132.12 133.33 134.55 135.76 136.97 138.18 139.39 140.61 141.82 143.03 144.24 145.45 146.67 147.88 149.09 150.3 151.52 152.73 153.94 155.15 156.36 157.58 158.79 160];
    angg = [134.5000    79.0000   58.1300   51.2500   51.2500   55.2900   62.1400   72.1400   86.7000  110.2800  120.0000  122.1100  124.2100  126.3200  128.4200  130.5300  132.6300  134.7400  136.8400  138.9500  141.0500  143.1600  145.2600  147.3700  149.4700  151.5800  153.6800  155.7900  157.8900  160.0000];
%     angg = [180	125	64.86	42.77	33.72	24.77	24.92	24.26	23.8	23.61	23.62	23.81	24.11	24.52	24.99	25.52	26.04	26.58	27.08	27.56	27.99	28.42	28.82	29.23	29.65	30.11	30.59	31.11	31.66	32.23	32.79	33.35	33.89	34.44	34.98	35.54	36.11	36.73	37.36	38.04	38.71	39.4	40.08	40.77	41.46	42.19	42.94	43.74	44.56	45.41	46.26	47.13	48.01	48.93	49.89	50.91	51.96	53.04	54.14	55.28	56.46	57.71	59.03	60.4	61.82	63.29	64.84	66.5	68.21	70.11	71.98	74.04	76.26	78.59	81.08	83.85	86.72	89.91	93.48	97.32	100];
    plot(1:length(angg),angg)
    angles = angg*pi/180;

    % ===== Simulate the Decay ======

    TR = ETL*TE+TRec;  % time until next excitation
    TEs= [];

    M = zeros(3,TR*nET*100,size(df,2));	% Keep track of magnetization at all time points. Acuracy 0.01ms

%     [Mprep]=T2Prep(T1,T2,T2ETL,T2int*TErep,M,df,b1);
% 
%     for ini=1:size(M,3)
%         M(:,1,ini)=Mprep(:,:,ini);	% Starting magnetization.
%     end

    for ini=1:size(M,3)
        M(:,1,ini)=[0;0;1];	% Starting magnetization.
    end

    % [M]=T2Prep(T1,T2,T2ETL,T2int,M,df,b1);

    for sp=1:size(M,3)
        for i=0:nET-1

            if(i==0) %excitation pulse
                M(:,1,sp)=excrot*M(:,1,sp);
                M(:,1,sp) = squeeze(K(sp,:,:))*M(:,1,sp)+squeeze(L(sp,:,:))';
            else
                M(:,i*TR*100+1,sp) = spoil*M(:,i*TR*100,sp);
                M(:,i*TR*100+1,sp) = excrot*M(:,i*TR*100+1,sp);
                M(:,i*TR*100+1,sp) = squeeze(K(sp,:,:))*M(:,i*TR*100+1,sp)+squeeze(L(sp,:,:))';
            end

            for k=0.02:0.01:TR %k = time (0.01ms)
                
                if(abs(mod(mod(k+i*TR,TR),TE)-TE/2)<0.000001 && k<ETL*TE) %RF pulse
                    plot(squeeze(M(1,:,:)+i*M(2,:,:)));
                    TEs = [TEs k+i*TR+TE/2];
                    cuberot = xrot(angles(round(mod((k+i*TR+TE/2),TR)/TE)));
                    M(:,round(k*100)+i*TR*100,sp) = squeeze(K(sp,:,:))*M(:,round(k*100)+i*TR*100-1,sp)+squeeze(L(sp,:,:))';
                    M(:,round(k*100)+i*TR*100,sp) = cuberot*M(:,round(k*100)+i*TR*100,sp); 
                elseif(mod(mod(k+i*TR,TR),TE)<TE/2 && k<ETL*TE) %free precession before RF
                    M(:,round(k*100)+i*TR*100,sp) = squeeze(K(sp,:,:))*M(:,round(k*100)+i*TR*100-1,sp)+squeeze(L(sp,:,:))';
                
                elseif((mod(mod(k+i*TR,TR),TE)>TE/2 && k<=ETL*TE) || k==ETL*TE) %free precession after RF
                    M(:,round(k*100)+i*TR*100,sp) = squeeze(K(sp,:,:))*M(:,round(k*100)+i*TR*100-1,sp)+squeeze(L(sp,:,:))';
                elseif(k>=ETL*TE) %free precession after echo train (recovery time)
                    M(:,round(k*100)+i*TR*100,sp) = squeeze(K(sp,:,:))*M(:,round(k*100)+i*TR*100-1,sp)+squeeze(L(sp,:,:))';

                end


            end
        end
    end


    % ===== Plot the Results ======
    clear i;

    %Calculate Complex Sum
    Mcomb=squeeze(M(1,:,:)+i*M(2,:,:)); 

    % figure('position',[200 200 1000 700])
    time = [0.01:0.01:nET*TR];
    % subplot(3,1,1)

    %plot Mx, My, Mz for each off resonance
    %plot phase for each off resonance
    % for pl=1:size(M,3)
    %     subplot(3,1,1)
    %     % plot magnitude:
    %     plot(time,(M(1,:,pl)),'b-',time,(M(2,:,pl)),'r-',time,M(3,:,pl),'k-','LineWidth',2);
    %     legend('M_x','M_y','M_z');
    %     hold on
    %     subplot(3,1,3)
    %     %plot phase:
    %     plot(time,angle(Mcomb(:,pl))) 
    %     hold on
    %     legend('phase')
    % end

    %Plot mean Mxy
    % subplot(3,1,2)
    if TErep == 1
        figure('position',[200 200 900 600]);
    end
    meanMxyMag = abs(mean(Mcomb,2));
    uTEs = unique(TEs);
    plot(time,abs(mean(Mcomb,2)),'b','LineWidth',2)
    plot(uTEs,meanMxyMag(round(uTEs*100)),'LineWidth',4)
    hold on
    scatter(TEs(1:nET*ETL),zeros(1,nET*ETL),'.k')

    legend('Mxy','Echo times');
    xlabel('Time (ms)');
    ylabel('Magnetization');
    axis([min(time) max(time) -0 0.6]);
    grid on;
    
     
    figure;
    [AX, H1, H2]=plotyy(uTEs,meanMxyMag(round(uTEs*100)),[uTEs], angg);
    set(get(AX(1),'Ylabel'),'String','MTF');
    set(get(AX(2),'Ylabel'),'String','Angle');
    set(AX(1),'YLim',[0 1])
    set(AX(1),'YTick',[0:0.1:1])
    set(AX(2),'YLim',[0 180]);
    set(AX(2),'YTick',[0:10:180])
    set(H1, 'LineWidth',4);
    set(H2, 'LineWidth',2,'Marker','o' );
    xlabel('Time (ms)');
    grid on;
    % title(['Flip angles x ',num2str(b1),', Number of Prepulses = ' num2str(T2ETL),', interval =  ',num2str(T2int),'ms'])
    meanmags(TErep)=mean(meanMxyMag(uTEs(5:end)));


