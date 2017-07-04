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

clear all

for TErep=1:6

    % ===== Set the parameters ======

    % Tissue properties
    nsp = 600;         % number of spins with different random off-resonance freqs
    T1 = 300;        % ms.
    T2 = 40;        % ms.
    df0 = 1;         % df0=1 will only display the results for one on resonance spin


    % Sequence Properties
    TE = 4;          % ms, repetition time of FSPGR
    ETL = 64;        % number of inversions
    TRec = 50;      % time to wait after echotrains
    nET = 1;         % number of echo trains
    invfa = 180;     % degree, flip angle inversion
    spoilpart = 0;   % [0 1],  0=perfect spoiling, 1=no spoiling, spoiling after echo train before next inversion.

    % T2 prep properties // 90x - 180y 180y (...) - -90x
    T2ETL = 4;      % number of refocussing pulses
    T2int = 6;      % time between refocussing pulses
    b1 =  180/180;       % factor multiplied with refocussing pulses (1 = 180 degree)


    % ===== Get the Propagation Matrices ======


    %calculate the off resonance frequencies [Hz]
    df = 1300*(rand(1,nsp)-0.5);  

    if df0
        df=0;
    end

    if size(df,2)==1
        df=[df df];
    end

    for fr=1:size(df,2) %free precession matrices for different dfs
        [K(fr,:,:),L(fr,:,:)] = freeprecess(1,T1,T2,df(fr));
    end

    invrot = xrot(invfa*2*pi/360);
    excrot = yrot(90*2*pi/360);


    spoil = xyspoil(spoilpart);

    % Get the flipangles for CUBE
    angles=CUBEangles(T1,T2,TE,0.24,ETL);
%     angles=CUBEangles(1200,80,TE,0.2,ETL);
    angles = b1 *angles;

    angdeg=angles*180/pi;
    % plot(1:length(angdeg),angdeg)

    % angg=[160 140 145 150 155 160:-4:20 20:4:160]
    % angg=[160 150 130 110 130 145 150 155 160:-4:30 28 27 25 24 22 21 20 21 22 22 23 24 25 26 28 30:4:120 126:6:160]
    % angg = [160 147 156 164 168 173 177 180*ones(1,121)];
    % plot(1:length(angg),angg)
    % angles = angg*pi/180;

    % ===== Simulate the Decay ======

    TR = ETL*TE+TRec;  % time until next excitation
    TEs= [];

    M = zeros(3,TR*nET,size(df,2));	% Keep track of magnetization at all time points.

    [Mprep]=T2Prep(T1,T2,T2ETL,T2int*TErep,M,df,b1);

    for ini=1:size(M,3)
        M(:,1,ini)=Mprep(:,:,ini);	% Starting magnetization.
    end

    % [M]=T2Prep(T1,T2,T2ETL,T2int,M,df,b1);

    for sp=1:size(M,3)
        for i=0:nET-1

            if(i==0) %excitation pulse
                M(:,1,sp)=excrot*M(:,1,sp);
                M(:,1,sp) = squeeze(K(sp,:,:))*M(:,1,sp)+squeeze(L(sp,:,:))';
            else
                M(:,i*TR+1,sp) = spoil*M(:,i*TR,sp);
                M(:,i*TR+1,sp) = excrot*M(:,i*TR+1,sp);
                M(:,i*TR+1,sp) = squeeze(K(sp,:,:))*M(:,i*TR+1,sp)+squeeze(L(sp,:,:))';
            end

            for k=2:TR %k = time (ms)
                if(mod(mod(k+i*TR,TR),TE)<TE/2 && k<ETL*TE) %free precession before RF
                    M(:,k+i*TR,sp) = squeeze(K(sp,:,:))*M(:,k+i*TR-1,sp)+squeeze(L(sp,:,:))';
                elseif(mod(mod(k+i*TR,TR),TE)==TE/2 && k<ETL*TE) %RF pulse
                    TEs = [TEs k+i*TR+TE/2];
                    cuberot = xrot(angles(mod((k+i*TR+TE/2),TR)/TE));
                    M(:,k+i*TR,sp) = squeeze(K(sp,:,:))*M(:,k+i*TR-1,sp)+squeeze(L(sp,:,:))';
                    M(:,k+i*TR,sp) = cuberot*M(:,k+i*TR,sp); 
                elseif((mod(mod(k+i*TR,TR),TE)>TE/2 && k<=ETL*TE) || k==ETL*TE) %free precession after RF
                    M(:,k+i*TR,sp) = squeeze(K(sp,:,:))*M(:,k+i*TR-1,sp)+squeeze(L(sp,:,:))';
                elseif(k>=ETL*TE) %free precession after echo train (recovery time)
                    M(:,k+i*TR,sp) = squeeze(K(sp,:,:))*M(:,k+i*TR-1,sp)+squeeze(L(sp,:,:))';

                end


            end
        end
    end


    %% ===== Plot the Results ======
    clear i;

    %Calculate Complex Sum
    Mcomb=squeeze(M(1,:,:)+i*M(2,:,:)); 

    % figure('position',[200 200 1000 700])
    time = [1:nET*TR];
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
    subplot(2,1,1)
    meanMxyMag = abs(mean(Mcomb,2));
    uTEs = unique(TEs);
     plot(time,abs(mean(Mcomb,2)),'b','LineWidth',2)
    plot(uTEs,meanMxyMag(uTEs),'LineWidth',4)
    hold on
     scatter(TEs(1:nET*ETL),zeros(1,nET*ETL),'.k')

    legend('Mxy','Echo times');
    xlabel('Time (ms)');
    ylabel('Magnetization');
    axis([min(time) max(time) -0 0.6]);
    grid on;
    % title(['Flip angles x ',num2str(b1),', Number of Prepulses = ' num2str(T2ETL),', interval =  ',num2str(T2int),'ms'])
    meanmags(TErep)=mean(meanMxyMag(uTEs(5:end)));

end

%Estimation
% estloc = 40;

% yd=meanmags';
% TE=[T2int*T2ETL]';
% for c=2:TErep
%     TE=[TE T2int*c*T2ETL];
% end
% 
% TE=TE';
% 
% initval = [10, 1, 0.01]';
% initR2s=5:3:25;
% initAs=[0:0.1:1.2];
% initvalmat=[10 0.1 0.01];
% 
% for rr=1:size(initR2s,2)
%     for aa=1:size(initAs,2)
%         initvalmat = [initvalmat; initR2s(rr) initAs(aa) 0.01];
%     end
% end
% initvalmat=initvalmat';
% 
% extraInitialValues = mat2cell(initvalmat,3,ones(1,size(initvalmat,2)));
% 
%     g = [ones(1,size(TE,1));-TE']';
%     initv =  g\log(yd);
%     initvalls = [initv(2)*1000 min(exp(initv(1)),1E6) 0.01]';
% 
% % Estimation
% fun = @(tht) predict_SET2( tht , TE);
% tht = fit_MRI( fun, yd, initvalls,'numPDFoptpar', 1); 
% % [CRLB, I, J] = CramerRaoLowerBound_MRI( tht(1:2,:,:), fun, tht(3,:,:));
% 
% 
% A = tht(2);
% T2e = 1000/tht(1);
% R2 = tht(1)/1000;
% 
% TEplot = [0 TE'];
% subplot(2,1,2)
% hold on
% % subplot(2,1,1)
% % scatter(estloc.*ones(size(TE)),yd,'.k','SizeData',450)
% % subplot(2,1,2)
% plot(TEplot,A.*exp(-TEplot*R2),'LineWidth',3)
% hold on
% scatter(TE,yd,'.k','SizeData',450)
% 
% 
% xlabel('Prepulse EchoTime (ms)')
% ylabel('Signal Intensity (au)')
% title('Fit through middle echo intensities')
% legend(['Estimated T2 = ' num2str(T2e) ' ms'])
% 
% 
% axis([0 max(TE)+10 0 0.35])



% 
% figure('position',[200 200 400 300]);
% ft=fft(meanMxyMag(uTEs));
% ft=ft';
% ftn=fliplr(ft);
% uTEsn = -fliplr(uTEs);
% fftt=[ftn ft];
% uTEsf=[uTEsn uTEs];
% plot(uTEsf,fftt)
% axis([-max(time)+300 max(time)-300 -3 0.3*max(fftt)]);
% 
% title(['Flip angles x ',num2str(b1),', Number of Prepulses = ' num2str(T2ETL),', interval =  ',num2str(T2int),'ms'])

% other plots/experiments:

% transmag=abs(mean(Mcomb,2));
% plot(TEs(1:nET*ETL),transmag(TEs(1:nET*ETL)))

% Phase0compMagn = abs(Mcomb(:,:)).*cos(angle(Mcomb(:,:)));
% plot(time,mean(Phase0compMagn,2),'k')


% ===== Do the estimations =====

% Set fit parameters.

% fiteq = ('a-b*exp(-x*c)');
% s = fitoptions('Method','NonLinearLeastSquares','Lower',[0,0,0],'Upper',[30000,30000,30000],'Startpoint',[1, 2, 0.01],'Maxiter',500,'Display','off','TolFun',10^-10,'TolX',10^-10);
% f = fittype(fiteq,'options',s);
% yd = M(3,(fitCurve-1)*nTR+1:(fitCurve)*nTR)';
% TI = time(1:nTR)';

% if TR can be varied flip angle can be estimated with this model
% fiteq = (['a-b*exp(-x*(c-log(cos(d*' num2str(pi/180) '))/' num2str(TR) '))'   ])
% to include TR between inversions, not necessary for correct T1 estimation:
% fiteq = (['(a-b*exp(-x*c)+exp(-' num2str(TRinv) '*c))' ]);

% Do the fitting
% [c2,gof,output] = fit(TI,yd,f,s);
% A=c2.a ;
% B=c2.b;
% R1=c2.c;
%     
% % Plot results 
% yy=A-B*exp(-TI*R1);
% subplot(2,1,2)
% plot(TI,yy,'g','LineWidth',3);
% hold on
% plot(TI,yd,':k','LineWidth',2);
% legend('Fitting result','Simulated Data','Location','SouthEast');
% 
% % Print out results
% T1s=1/R1;
% T1cor = (1/T1s + (1/TR)*log(cos(fa*2*pi/360)))^-1;
% T1cor2 = T1s * (B/A - 1);
% disp(['T1* = ' num2str(T1s) ' ms'])
% disp(['corrected T1 (flip angle/TR correction)= ' num2str(T1cor) ' ms'])
% disp(['2nd corrected T1 (B/A correction) = ' num2str(T1cor2) ' ms'])
% 
% % disp(['A = ' num2str(A) ' '])
% % disp(['B = ' num2str(B) ' '])
% 
% % disp(['a = ' num2str(a) ' '])
% 
% % ML and Cramer-Rao estimations for later implementation
% % yd = abs(M(3,3*nTR+1:4*nTR))';
% % TI = time(3*nTR+1:4*nTR)';
% % CRLBsigma = 1;
% % opt = optimset('fminunc');
% % opt = optimset(opt,'Diagnostics','off','LargeScale','off','gradObj','on','Display','off','MaxIter',80,'Hessian','off','TolFun',1e-12,'Tolx',1e-8);
% % s = fitoptions('Method','NonLinearLeastSquares','Lower',[0,0,0],'Upper',[30000,30000,100],'Startpoint',[8200, 0.001, 0.001],'Maxiter',500,'Display','off','TolFun',10^-5,'TolX',10^-5);
% % fiteq = (['a*(1-b*exp(-x*c)+ exp(-(' num2str(TR(1)) '+x)*c))' ]);
% % f = fittype(fiteq,'options',s);
% % fitrange = 5:size(TI);
% % [LS,ML] = T1IRAbsComputePar(yd,TI,0,CRLBsigma,opt,s,f, TR,fitrange);
% % [CR] = T1IRAbsCramerRao(ML,TI,CRLBsigma, TR);
