% Bloch Equation Simulation for mCINE-IR
%
% Thanks to Brian Hargreaves, Stanford,
% http://mrsrl.stanford.edu/~brian/bloch/
%
% Henk Smit, h.smit@erasmusmc.nl, 03-2012.
% -----------------------------------------

clear all

bp = 375:25:375;
bp = repmat(bp,1,1);
RRs = round(60000./bp);

for rr=1:size(RRs,2)

    % ===== Set the parameters ======

    % Tissue 1 properties
    df = 0;         % Hz off-resonance
    T1 = 900:200:1100;      % ms.
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
    sigma = 0.00;

    % Sequence Properties
    TR = 15;     % ms, repetition time of FSPGR
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
%     M=abs(M+noise);
    M=abs(M);
    
    noisem=sigma*(randn(size(ydmix,1),size(ydmix,2))+i*randn(size(ydmix,1),size(ydmix,2)));
    ydmix=abs(ydmix+noisem);
%     ydmix=abs(ydmix);
    
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


    yd = squeeze(M(3,1,(fitCurve-1)*nTR+1:fitCurve*nTR));
    yd2 = squeeze(M(3,2,(fitCurve-1)*nTR+1:fitCurve*nTR));
    TI = time(1:nTR)';



    % if TR can be varied flip angle can be estimated with this model
    % fiteq = (['a-b*exp(-x*(c-log(cos(d*' num2str(pi/180) '))/' num2str(TR) '))'   ])
    % to include TR between inversions, not necessary for correct T1 estimation:
    % fiteq = (['(a-b*exp(-x*c)+exp(-' num2str(TRinv) '*c))' ]);

    % Do the fitting
%     [c,gof,output] = fit(TI,yd,f,s);
%     A=c.a ;
%     B=c.b;
%     R1=c.c;
% 
%     [c2,gof,output] = fit(TI,yd2,f,s);
%     A2=c2.a ;
%     B2=c2.b;
%     R12=c2.c;
% 
%     [cm,gof,output] = fit(TI,ydmix,f,s);
%     Am=cm.a ;
%     Bm=cm.b;
%     R1m=cm.c;

initval = [2, 1, 2, 0.1]';
initR1s=0:0.15:6;
initAs=0.7:0.1:1;
initBs=1.7:0.1:2;
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


% Estimation
fun = @(tht) predict_IRT1( tht , TI);
tht = fit_MRI( fun, yd, initval,'numPDFoptpar', 1,'initialValueSpecifierVect',[1 0 0 1],'extraInitialValues',extraInitialValues); 
[CRLB, I, J] = CramerRaoLowerBound_MRI( tht(1:3,:,:), fun, tht(4,:,:));

ML(1) = tht(2);
ML(2) = tht(3);
ML(3) = tht(1);
ML=ML';
    
A = ML(1);
B = ML(2);
T1 = 1000/ML(3);
R1 = ML(3)/1000;

% fun = @(tht) predict_IRT1( tht , TI);
% tht2 = fit_MRI( fun, yd2, initval,'numPDFoptpar', 1,'initialValueSpecifierVect',[1 0 0 1],'extraInitialValues',extraInitialValues,'noiseLevel',0.001); 
% [CRLB2, I2, J2] = CramerRaoLowerBound_MRI( tht(1:3,:,:), fun, tht(4,:,:));

initval2 = [2, 1, 2]';
fun = @(tht) predict_IRT1( tht , TI);
tht2 = fit_MRI( fun, yd2, initval2,'numPDFoptpar', 0,'noiseLevel',0.000001); 
[CRLB2, I2, J2] = CramerRaoLowerBound_MRI( tht(1:3,:,:), fun, tht(4,:,:));

ML2(1) = tht2(2);
ML2(2) = tht2(3);
ML2(3) = tht2(1);
ML2=ML2';
    
A2 = ML2(1);
B2 = ML2(2);
T12 = 1000/ML2(3);
R12 = ML2(3)/1000;

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
    yy=A2-B2*exp(-TI*R12);
%     subplot(3,1,2)
    plot(TI,yy,'g','LineWidth',3);
    hold on
    plot(TI,yd2,':k','LineWidth',2);
    legend('Fitting result','Simulated Data','Location','SouthEast');
% 
%     plot(TI,yd,'r--','LineWidth',1)
%     hold on
%     plot(TI,yd2,'r--','LineWidth',1)
%     plot(TI,ydmix,'k','LineWidth',2)

    % Print out results
    T1est(rr)=1/R1;
    T1est2(rr)=1/R12;
%     T1mest(rr)=1/R1m;

end

bpms = 60000./RRs;
plot(bpms,T1est,'--g','LineWidth',3)
hold on
plot(bpms,T1est2,'--b','LineWidth',3)
plot(bpms,T1mest,':r','LineWidth',3)

scatter(bpms,T1est,'*g')
hold on
scatter(bpms,T1est2,'*b')

xlabel('Heart Rate (bpm)')
ylabel('Estimated T1 (ms')
legend('Tissue 1','Tissue 2','Mixed tissue')
axis([min(bpms)-15 max(bpms)+15 0 1450]);

plot([100 800],[242 242],'k','LineWidth',1)
plot([100 800],[1043 1043],'k','LineWidth',1)


yy=abs(A-B.*exp(-R1.*TI));
plot(TI,yy,'r')
hold on
scatter(TI,yd,'.r')

yy2=abs(A2-B2.*exp(-R12.*TI));
plot(TI,yy2,'b')
hold on
scatter(TI,yd2,'.b')

yym=abs(Am-Bm.*exp(-R1m.*TI));
plot(TI,yym,'k')
hold on
scatter(TI,ydmix,'.k')


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
