function [Mprep] = T2prep(T1,T2,ETL,TE,Ms,df,b1)
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



% ===== Set the parameters ======

% Tissue properties
nsp = size(Ms,3);         % number of spins with different random off-resonance freqs


% Sequence Properties
TRec = 1;      % time to wait after echotrains
% nET = 1;         % number of echo trains
% invfa = 180;     % degree, flip angle inversion
% spoilpart = 0;   % [0 1],  0=perfect spoiling, 1=no spoiling, spoiling after echo train before next inversion.

% ===== Get the Propagation Matrices ======


%calculate the off resonance frequencies [Hz]

for fr=1:size(df,2) %free precession matrices for different dfs
    [K(fr,:,:),L(fr,:,:)] = freeprecess(1,T1,T2,df(fr));
end

excrot = yrot(90*2*pi/360);
excrot2 = yrot(-90*2*pi/360);

% spoil = xyspoil(spoilpart);

% Set the refocussing flip angles 
angles=pi*ones(1,ETL);
angles = b1*angles;

% ===== Simulate the Decay ======

TR = ETL*TE+TRec;  % time until next excitation
TEs= [];
nET = 1;

M = zeros(3,TR*nET,size(df,2));	% Keep track of magnetization at all time points.

for ini=1:size(M,3)
	M(:,1,ini)=[0;0;1];	% Starting magnetization.
end


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
            elseif(mod(mod(k+i*TR,TR),TE)==TE/2 && k<ETL*TE) %RF pulse if mod k = TE/2
                TEs = [TEs k+i*TR+TE/2];
                cuberot = xrot(angles(mod((k+i*TR+TE/2),TR)/TE));
                M(:,k+i*TR,sp) = squeeze(K(sp,:,:))*M(:,k+i*TR-1,sp)+squeeze(L(sp,:,:))';
                M(:,k+i*TR,sp) = cuberot*M(:,k+i*TR,sp); 
            elseif((mod(mod(k+i*TR,TR),TE)>TE/2 && k<=ETL*TE) || k==ETL*TE) %free precession after RF
                M(:,k+i*TR,sp) = squeeze(K(sp,:,:))*M(:,k+i*TR-1,sp)+squeeze(L(sp,:,:))';
            elseif(k>=ETL*TE) %free precession after echo train (recovery time)
                M(:,k+i*TR,sp) = excrot2*M(:,k+i*TR-1,sp);

            end


        end
    end
end

Mprep = M(:,TR,:);


% ===== Plot the Results ======
% clear i;
% 
% %Calculate Complex Sum
% Mcomb=squeeze(M(1,:,:)+i*M(2,:,:)); 
% 
% figure('position',[200 200 1000 700])
% time = [1:nET*TR];
% subplot(3,1,1)
% hold on
% 
% %plot Mx, My, Mz for each off resonance
% %plot phase for each off resonance
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
% 
% %Plot mean Mxy
% subplot(3,1,2)
% meanMxyMag = abs(mean(Mcomb,2));
% uTEs = unique(TEs);
% % plot(time,abs(mean(Mcomb,2)),'b','LineWidth',2)
% plot(uTEs,meanMxyMag(uTEs))
% hold on
% scatter(TEs(1:nET*ETL),zeros(1,nET*ETL),'.k')
% 
% legend('Mxy','Echo times');
% xlabel('Time (ms)');
% ylabel('Magnetization');
% axis([min(time) max(time) -1.1 1.1]);
% grid on;

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
