% Bloch Equation Simulation for mCINE-IR
%
% Thanks to Brian Hargreaves, Stanford,
% http://mrsrl.stanford.edu/~brian/bloch/
%
% Henk Smit, h.smit@erasmusmc.nl, 03-2012.
% -----------------------------------------

clear all

bp = 365:20:365;
bp = repmat(bp,1,5);
RRs = round(60000./bp);

for rr=1:size(RRs,2)
    
      disp(['Calculating : ' num2str(rr) '/' num2str(size(RRs,2))]);

    % ===== Set the parameters ======

    % Tissue 1 properties
    df = 0;         % Hz off-resonance
    T1 = [450 1150 2400];      % ms.
    T2 = 50;       % ms.

    % Tissue 2 properties
    df = 0;         % Hz off-resonance
    T12 = 400;      % ms.
    T22 = 50;       % ms.
    
    % Tissue ratios:
    T1rat = 1;
    T2rat = 1;

    % Heartrate properties. Heartrate changes to be implemented!
    nRR = 15;       % number of RR intervals between inversions
    RR = RRs(rr);  % ms, RR interval
    
    % Noise sigma
    sigma = 0.03;

    % Sequence Properties
    TR = 3.8;     % ms, repetition time of FSPGR
    fa = 3.5;         % degree, flip angle FSPGR
    nTI = 3;        % number of inversions
    invfa = 180;    % degree, flip angle inversion
    spoilpart = 0;  % [0 1], part of transversal magnetization that is is not spoiled. 0=perfect spoiling, 1=no spoiling

    % Estimation properties
    fitCurve = 3;   % which curve is used for the estimation. For the first curve both correction factors work.
    
    



    % ===== Get the Propagation Matrices ======

    nT1 = size(T1,2);
    nBPM = size(RR,2);
    Kcell = cell(nT1,1);
    Lcell = cell(nT1,1);

    for tt=1:nT1
        [K,L] = freeprecess(TR,T1(tt),T2,df);
        Kcell(tt)={K};
        Lcell(tt)={L};
    end

    invrot = throt(invfa*2*pi/360,0);

    spgrrot = throt(fa*2*pi/360,pi);

    spoil = xyspoil(spoilpart);


    % ===== Simulate the Decay ======


    nTR = floor((nRR)*RR/TR);
    TRinv = nRR*RR; 
 

    M = zeros(3,nT1,nTR*nTI);	% Keep track of magnetization at all time points.
    for nt1=1:nT1
        M(:,nt1,1)=[0;0;1];	% Starting magnetization.
    end
    clear nt1;


    for nt1=1:nT1    %different tissues
        for j=0:nTI-1

            if(j==0)
                M(:,nt1,1)=invrot*M(:,nt1,1);
            else
                M(:,nt1,j*nTR+1) = invrot*M(:,nt1,(j*nTR));
            end

            for k=2:nTR

                M(:,nt1,j*nTR+k) = Kcell{nt1}*M(:,nt1,j*nTR+k-1)+Lcell{nt1};
                M(:,nt1,j*nTR+k) = spoil*M(:,nt1,j*nTR+k);
                M(:,nt1,j*nTR+k) = spgrrot*M(:,nt1,j*nTR+k);

            end
        end
    end
    
    ydmix=mean([T1rat.*squeeze(M(3,1,(fitCurve-1)*nTR+1)) T2rat.*squeeze(M(3,2,(fitCurve-1)*nTR+1)) ]);
    for k=(fitCurve-1)*nTR+2:fitCurve*nTR
        ydmix=[ydmix mean([T1rat.*M(3,1,k) T2rat.*M(3,2,k)])];
    end
    ydmix=ydmix';
    
    noise=sigma*(randn(size(M,1),size(M,2),size(M,3))+i*randn(size(M,1),size(M,2),size(M,3)));
    M=abs(M+noise);
    M=abs(M);
    
    noisem=sigma*(randn(size(ydmix,1),size(ydmix,2))+i*randn(size(ydmix,1),size(ydmix,2)));
    ydmix=abs(ydmix+noisem);
    ydmix=abs(ydmix);
    
    clear noise
    clear noisem
    
    

    % ===== Plot the Results ======

    time = [0:nTI*nTR-1]*TR;
%     subplot(3,1,1)
%     hold on
% 
%     for ip=1:nT1
%         plot(time,(squeeze(M(1,ip,:))),'b-',time,(squeeze(M(2,ip,:))),'r--',time,squeeze(M(3,ip,:)),'r-','LineWidth',2);
%         hold on
%     end
% 
%     legend('M_x','M_y','M_z');
%     xlabel('Time (ms)');
%     ylabel('Magnetization');
%     axis([min(time) max(time) -1 1]);
%     grid on;
%     plot(time,(M2(1,:)),'b-',time,(M2(2,:)),'r--',time,M2(3,:),'k-','LineWidth',2);


    % ===== Do the estimations =====

    % Set fit parameters.

%     fiteq = ('a-b*exp(-x*c)');
%     s = fitoptions('Method','NonLinearLeastSquares','Lower',[0,0,0],'Upper',[30000,30000,30000],'Startpoint',[1, 2, 0.01],'Maxiter',500,'Display','off','TolFun',10^-10,'TolX',10^-10);
%     f = fittype(fiteq,'options',s);

clear yd
for sigs=1:nT1
    yd(:,sigs) = squeeze(M(3,sigs,(fitCurve-1)*nTR+1:fitCurve*nTR));
    TI = time(1:nTR)';
end


initval = [2, 1, 2, 0.1]';
initR1s=1000./T1;
initAs=0.4:0.2:0.8;
initBs=0.8:0.2:1.7;
initvalmat=[2 1 2 0.1];

for initq=1:size(initR1s,2)
    for aa=1:size(initAs,2)
        for bb=1:size(initBs,2)
            initvalmat = [initvalmat; initR1s(initq) initAs(aa) initBs(bb) 0.1];
        end
    end
end
initvalmat=initvalmat';
extraInitialValues = mat2cell(initvalmat,4,ones(1,size(initvalmat,2)));
clear initvalmat

for dT1=1:nT1

% Estimation
    fun = @(tht) predict_IRT1( tht , TI);
    tht = fit_MRI( fun, yd(:,dT1), initval,'numPDFoptpar', 1,'initialValueSpecifierVect',[1 0 0 1],'extraInitialValues',extraInitialValues); 
    [CRLB, I, J] = CramerRaoLowerBound_MRI( tht(1:3,:,:), fun, tht(4,:,:));

    ML(1) = tht(2);
    ML(2) = tht(3);
    ML(3) = tht(1);
    ML=ML';

    A(dT1) = ML(1);
    B(dT1) = ML(2);
    T1(dT1) = 1000/ML(3);
    R1(dT1) = ML(3)/1000;
end


% fun = @(tht) predict_IRT1( tht , TI);
% thtm = fit_MRI( fun, ydmix, initval,'numPDFoptpar', 1,'initialValueSpecifierVect',[1 0 0 1],'extraInitialValues',extraInitialValues); 
% [CRLBm, Im, Jm] = CramerRaoLowerBound_MRI( tht(1:3,:,:), fun, tht(4,:,:));
% 
% MLm(1) = thtm(2);
% MLm(2) = thtm(3);
% MLm(3) = thtm(1);
% MLm=MLm';
%     
% Am = MLm(1);
% Bm = MLm(2);
% T1m = 1000/MLm(3);
% R1m = MLm(3)/1000;



    % Plot results 
%     yy=A2-B2*exp(-TI*R12);
% %     subplot(3,1,2)
%     plot(TI,yy,'g','LineWidth',3);
%     hold on
%     plot(TI,yd2,':k','LineWidth',2);
%     legend('Fitting result','Simulated Data','Location','SouthEast');
% 
%     plot(TI,yd,'r--','LineWidth',1)
%     hold on
%     plot(TI,yd2,'r--','LineWidth',1)
%     plot(TI,ydmix,'k','LineWidth',2)

    % Print out results
    T1est(rr,:)=1./R1;
    Aest(rr,:)=A;
    Best(rr,:)=B;

%     T1mest(rr)=1/R1m;

end

a = figure('position',[200 200 500 400]);

for ii=1:size(T1est,2)
    scatter(TI,yd(:,ii),'.','MarkerEdgeColor',(ii-1).*[0.2 0.2 0.2]+[0.3 0.3 0.3])
    hold on
end
legend('T1 = 369ms','T1 = 735ms','T1 = 1107ms','Location','NorthWest')
for ii=1:size(T1est,2)
    plot(TI,abs(mean(Aest(:,ii))-mean(Best(:,ii)).*exp(-TI/mean(T1est(:,ii)))),'--k','LineWidth',2)
end

xlabel('Inversion Time (ms)','FontSize', 13)
ylabel('Signal Intensity (au)','FontSize', 13)
% legend('Tissue 1','Tissue 2','Mixed tissue')
axis([0 2500 0 1.1]);

% set(gca,'fontsize',9)
% name = ['some_name9.tif'];
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3.1])
% print(a, '-r320', '-dtiff', name);

% yy=abs(A-B.*exp(-R1.*TI));
% plot(TI,yy,'r')
% hold on
% scatter(TI,yd,'.r')
% 
% yy2=abs(A2-B2.*exp(-R12.*TI));
% plot(TI,yy2,'b')
% hold on
% scatter(TI,yd2,'.b')
% 
% yym=abs(Am-Bm.*exp(-R1m.*TI));
% plot(TI,yym,'k')
% hold on
% scatter(TI,ydmix,'.k')


% CRLBsigma = (size(yd,1)-fitfrom)*tht(4)/(size(yd,1)-size(tht,1)-fitfrom+1);
% CR = zeros(3,3);
% CR(3,3)=CRLB(1); 
% CR(1,1)=CRLB(3);
% CR(2,2)=CRLB(6);




% 
% fiteqlin = ('a*x+c');
% slin = fitoptions('Method','NonLinearLeastSquares','Lower',[-10,-10],'Upper',[30000,30000],'Startpoint',[0, 0],'Maxiter',500,'Display','off','TolFun',10^-10,'TolX',10^-10);
% flin = fittype(fiteqlin,'options',slin);
% [clin,gof,output] = fit(bpms',T1mest',flin,slin);
% 
% alin=clin.a;
% blin=clin.c;
% range=200:600;
% plot(range,alin.*range+blin);
    

    



% ML and Cramer-Rao estimations for later implementation
% yd = abs(M(3,3*nTR+1:4*nTR))';
% TI = time(3*nTR+1:4*nTR)';
% CRLBsigma = 1;
% opt = optimset('fminunc');
% opt = optimset(opt,'Diagnostics','off','LargeScale','off','gradObj','on','Display','off','MaxIter',80,'Hessian','off','TolFun',1e-12,'Tolx',1e-8);
% s = fitoptions('Method','NonLinearLeastSquares','Lower',[0,0,0],'Upper',[30000,30000,100],'Startpoint',[8200, 0.001, 0.001],'Maxiter',500,'Display','off','TolFun',10^-5,'TolX',10^-5);
% fiteq = (['a*(1-b*exp(-x*c)+ exp(-(' num2str(TR(1)) '+x)*c))' ]);
% f = fittype(fiteq,'options',s);
% fitrange = 5:size(TI);
% [LS,ML] = T1IRAbsComputePar(yd,TI,0,CRLBsigma,opt,s,f, TR,fitrange);
% [CR] = T1IRAbsCramerRao(ML,TI,CRLBsigma, TR);
