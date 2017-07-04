%function [A,CRA,R2,CRR2] =simestimate(image, TE, sigma)

%simulate an A and R2 estimating experiment
%Henk Smit
%12.01.2011

clear('all')

%the 5 parameters you can change for a simulation, with in [] the value they
%have now

TE=0.001.*[1.3 2 3 5 7.3]'; %the echotimes in ms [5,15,..55]
R2real=[400]; %the R2's of the simdata [10,20...70]
As=1250;%the A's of the simdata [200,300..900]
R2realmore=R2real(ones(1,200),:); %the size of each A/R2 cluster(number of repetitions) [50]
sigma=20; %the sigma of the noise in the data [9]
Cs=0;

R2real=reshape(R2realmore,1,numel(R2realmore));

nroccur=size(R2realmore,1);%
nrR2s=size(R2real,2)/nroccur;
nrTEs=size(TE,1);
nrAs=size(As,2);
nrpointsperR2=nroccur*nrR2s;
nrclusters=nrAs*nrR2s;
nrpoints=nrAs*nrR2s*nroccur;


Amat=As(ones(1,nrTEs),:);

for g=0:(size(As,2)-1)
    image(:,1+(g*nrpointsperR2):(g+1)*nrpointsperR2) = Amat(:,(g+1)*ones(1,nrpointsperR2)).*(exp(-TE*R2real))+Cs;
end

%Noise can be added by enabling the following three lines
noise=sigma*(randn(size(image,1),size(image,2))+i*randn(size(image,1),size(image,2)));
image=image+noise;
image=abs(image);

sizey=size(image,2);
sizex=size(image,1);

A=zeros(sizey,1);
CRA=A;
R2=A;
CRR2=A;

opt = optimset('fminunc');
opt = optimset(opt,'LargeScale','off','gradObj','on','DerivativeCheck','off','Display','off','MaxIter',200,'Hessian','off','TolFun',1e-10,'TolX',1e-7,'MaxFunEvals',1000);
s = fitoptions('Method','NonLinearLeastSquares','Lower',[-10000,-1000,0],'Upper',[50000,10000,10000],'Startpoint',[1000,0,0],'Maxiter',12,'TolFun',10^-9,'TolX',10^-9,'Display','off');
f = fittype('a*exp(-x*b)+c','options',s);
fitrange = 1:size(TE);

for j=1:sizey
    disp(['Calculating pixel: ' num2str(j) '/' num2str(sizey)]);
          for k=1:sizex
             yd(k)=double(image(k,j));
          end

            %Possibly LS is not the really the LEAST squares solution!
            %
            [LS,ML]=T2ABCComputePar(yd',TE,0,sigma,opt,s,f,fitrange);
            [CR]=T2ABCCramerRao(ML,TE,sigma);
            
            DT_LS(1,1)=LS(1);
            DT_LS(1,2)=LS(2);
            DT_LS(1,3)=LS(3);
            
            R2(j,1) = ML(2,1);
            A(j,1) = ML(1,1);
            C(j,1) = ML(3,1);
            T2(j,1) = 1000/ML(2,1);
            
            R2ls(j,1)=DT_LS(1,2);
            Als(j,1)=DT_LS(1,1);
            Cls(j,1)=DT_LS(1,3);
            T2ls(j,1) = 1000/DT_LS(1,2);
            
            CRA(j,1) = CR(1,1);
            CRR2(j,1) = CR(2,2);  
            STDCRA(j,1) = sqrt(CR(1,1));
            STDCRT2(j,1) = 1000*sqrt(CR(2,2))/(ML(2,1)^2);
            
            
 end  

snr=zeros(nrclusters,1);

for l=0:(nrclusters-1)   
    meanr2=mean(R2((l*nroccur)+1:(l+1)*nroccur));
    stdr2=std(R2((l*nroccur)+1:(l+1)*nroccur));
    snr(l+1,1)=meanr2/stdr2;
    %snr(l+1,2)=R2real(l*nroccur);
end

snrR2=R2realmore(1,:);
snrA=As';
snr=reshape(snr,nrR2s,nrAs);

std(T2);
1000*mean(STDCRT2);

 scatter(TE,image(:,1),'.b')
 hold on;
 te=0:0.001:(max(TE)+0.02);
%  yy=ML(1,1)*exp(-te*ML(2,1))+ML(3,1);
%  plot(te,yy,'k');
 

 scatter(reshape(repmat(TE,1,nrpoints),sizex*sizey,1),reshape(image,sizex*sizey,1),'.b')
 hold on;
 te = 0:0.001:(max(TE)+0.02);
 yy = median(A)*exp(-te*median(R2))+median(C);
 yyls = median(Als)*exp(-te*median(R2ls))+median(Cls);
 yytrue = As*exp(-te*snrR2)+Cs;
 plot(te,yy,'k');
 plot(te,yyls,'r');
 plot(te,yytrue,'y');
 
    disp(['Model = A * Exp(-R2*TE) + B'])
    disp(['Results: medians and true values'])
    disp(['True R2 = ' num2str(snrR2) '      T2 = ' num2str(1000/snrR2) '  A = ' num2str(As) '      C = ' num2str(Cs)])
    disp(['ML   R2 = ' num2str(median(R2)) ' T2 = ' num2str(1000/median(R2)) ' A = ' num2str(median(A)) ' C = ' num2str(median(C))])
    disp(['LS   R2 = ' num2str(median(R2ls)) '  T2 = ' num2str(1000/median(R2ls)) '  A = ' num2str(median(Als)) ' C = ' num2str(median(Cls))])
    
%To show an A vs. R2 map with errorbars(either parallel with x and y axis
%or in the direction of the eigenvectors of the covariancematrix) enable one of the following lines: 

%errorbarxy(nonzeros(A)',nonzeros(T2)',sqrt(nonzeros(CRA))',nonzeros(STDCRT2)',[],[],'w')
%errorbarxyoblique(nonzeros(A)',nonzeros(R2)',lx',ly',ux',uy','w')
% hold on
% scatter(nonzeros(A),nonzeros(T2),'.k')
% scatter([400 600 800 400 600 800 400 600 800],[30 30 30 70 70 70 110 110 110],'dk','SizeData',75,'MarkerFaceColor','yellow');
% 
% xlabel('A');
% ylabel('T2');

%enable this line to create a surface plot of the snr a function of A and R
%surf (snrA,snrR2,snr, 'DisplayName', 'snrA'); 
% 
% scatter(TE(1*ones(size(image,2),1)),image(1,:),'*b')
% hold on
% te=0:0.002:max(TE)+min(TE);
% siR2 = @(A,te,R2) A*exp(-te*R2);
% siR2abc = @(A,te,R2,C) A*exp(-te*R2)+C;
% 
% % yy=siR2(mean(A(:,1)),te,mean(R2(:,1)));
% % yyls=siR2(mean(Als2(:,1)),te,mean(R2ls2(:,1)));
% % yytrue=siR2(As,te,snrR2);
% 
% % plot(te,yy(1,:),'r','LineWidth',1);
% % hold on;
% % plot(te,yyls(1,:),'b','LineWidth',1);
% 
% % hold on
% % plot(te,yytrue,'-wo','LineWidth',2, 'MarkerEdgeColor','k','MarkerFaceColor',[0.99 0.98 0.98], 'MarkerSize',5);
% % scatter(te,yytrue,'o','MarkerEdgeColor','k','MarkerFaceColor',[0.99 0.98 0.98]);
% % plot(te,yy(1,:),'r','LineWidth',1);
% 
% yyabc=siR2abc(mean(A(:,1)),te,mean(R2(:,1)),mean(C(:,1)));
% yylsabc=siR2abc(mean(Als2(:,1)),te,mean(R2ls2(:,1)),mean(Cls2(:,1)));
% yytrueabc=siR2abc(As,te,snrR2,As./As);
% 
% plot(te,yyabc(1,:),'r','LineWidth',1);
% hold on;
% plot(te,yylsabc(1,:),'b','LineWidth',1);
% 
% hold on
% plot(te,yytrueabc,'-wo','LineWidth',2, 'MarkerEdgeColor','k','MarkerFaceColor',[0.99 0.98 0.98], 'MarkerSize',5);
% scatter(te,yytrueabc,'o','MarkerEdgeColor','k','MarkerFaceColor',[0.99 0.98 0.98]);
% plot(te,yyabc(1,:),'r','LineWidth',1);
% 
% legend('Simulated Data','ML Fit','LS Fit','True Decay');
% 
% for i=2:size(image,1)
%    scatter(TE(i*ones(size(image,2),1)),image(i,:),'*b')
%    hold on;
% end
% 
% xlabel('time (s)');
% ylabel('signal intensity');



% plot(te,yytrueabc,'-wo','LineWidth',2, 'MarkerEdgeColor','k','MarkerFaceColor',[0.99 0.98 0.98], 'MarkerSize',5);
% hold on
% plot(te,yyabc,'--k','LineWidth',2);

%t2mat=[T2 T2ls]
%plot([0 100],[14.2 14.2],':k')
%hold on
%boxplot(t2mat,'labels',['ML Estimation';'LS Estimation'])
%legend('True T2 Value','Location','NorthWest');



%LSSI=siR2(DT_LS(1),TE,DT_LS(2));
%data=image(:,size(image,2));
%LSSI-data;
%ans.*ans;
%sum(ans);

%MLSI=siR2(ML(1),TE,ML(2));
%data=image(:,size(image,2));
%MLSI-data;
%ans.*ans;
%sum(ans);

