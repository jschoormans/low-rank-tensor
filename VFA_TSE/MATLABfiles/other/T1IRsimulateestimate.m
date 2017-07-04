%function [A,CRA,R1,CRR1] =simestimate(image, TE, sigma)

%simulate an A and R1 estimating experiment
%Henk Smit
%12.01.2011

clear('all')

%the 5 parameters you can change for a simulation, with in [] the value they
%have now

TI=([100:350:2200])';%([100 200 400 800 2100])'; %the inversiontimes in ms [5,15,..55]
TR=4*ones(1,size(TI,1))';%([200 200 200 200 200])';
R1real=0.0007; %the R1's of the simdata [10,20...70]
Bs=2;%the Bs's of the simdata [200,300..900]
As=7700; %The A of the simdata
R1realmore=R1real(ones(1,20),:); %the size of each A/R1 cluster(number of repetitions) [50]
sigma=100; %the sigma of the noise in the data [9]

R1real=reshape(R1realmore,1,numel(R1realmore));

nroccur=size(R1realmore,1);%
nrR1s=size(R1real,2)/nroccur;
nrTEs=size(TI,1);
nrBs=size(Bs,2);
nrpointsperR1=nroccur*nrR1s;
nrclusters=nrBs*nrR1s;
nrpoints=nrBs*nrR1s*nroccur;


Bmat=Bs(ones(1,nrTEs),:);


for g=0:(size(Bs,2)-1)
%     image(:,1+(g*nrpointsperR1):(g+1)*nrpointsperR1)=As-Bmat(:,(g+1)*ones(1,nrpointsperR1)).*(exp(-TI*R1real));
%     image2(:,1+(g*nrpointsperR1):(g+1)*nrpointsperR1)=(0.3*As)-(0.3*Bmat(:,(g+1)*ones(1,nrpointsperR1))).*(exp(-TI*1.25*R1real));
    image(:,1+(g*nrpointsperR1):(g+1)*nrpointsperR1)=As*(1-Bmat(:,(g+1)*ones(1,nrpointsperR1)).*(exp(-TI*R1real))+exp(-(TR+TI).*1*R1real));
end

% image=image+image2; %biexponential

%Noise can be added by enabling the following three lines
noise=sigma*(randn(size(image,1),size(image,2))+i*randn(size(image,1),size(image,2)));
image=abs(image+noise);
image=abs(image);

sizey=size(image,2);
sizex=size(image,1);

A=zeros(sizey,1);
CRA=A;
R1=A;
R1ls=A;
CRR1=A;
B=A;
Als=A;
Bls=A;

opt = optimset('fminunc');
opt = optimset(opt,'Diagnostics','off','LargeScale','off','DerivativeCheck','off','gradObj','on','Display','off','MaxIter',200,'Hessian','off','TolFun',1e-10,'Tolx',1e-8,'MaxFunEvals',200);%'DiffMaxChange',0.01);
s = fitoptions('Method','NonLinearLeastSquares','Lower',[0,0,0],'Upper',[30000,30000,100],'Startpoint',[3000, 2, 0.001],'Maxiter',10,'Display','off','TolFun',10^-5,'TolX',10^-5);
% f = fittype('a-(b*exp(-x*c))','options',s);
fiteq = (['a*(1-b*exp(-x*c)+ exp(-(' num2str(TR(1)) '+x)*c))' ])
% f = fittype('a*(1-b*exp(-x*c)+0.1422)','options',s);
f = fittype(fiteq,'options',s);
fitrange=2:size(TI);


for j=1:sizey
    disp(['Calculating pixel: ' num2str(j) '/' num2str(sizey)]);
          for k=1:sizex
             yd(k)=double(image(k,j));
          end
            [LS,ML]=T1IRAbsComputePar(yd',TI,0,sigma,opt,s,f,TR,fitrange);
            [CR]=T1IRAbsCramerRao(ML,TI,sigma,TR);

            %[LS,ML]=T1IRComputePar(yd',TI,[],0,sigma)
            %[CR]=T1IRCramerRao(ML,TI,sigma);
            %[EigVec(:,:,j),EigVal(:,:,j)]=eig(CR);
            
            R1(j,1) = ML(3,1);
            B(j,1) = ML(2,1);
            A(j,1) = ML(1,1);
           
            Als(j,1) = LS(1,1);
            R1ls(j,1) = LS(3,1);
            Bls(j,1) = LS(2,1);
            
            CRA(j,1) = CR(1,1);
            CRB(j,1) = CR(2,2);
            CRR1(j,1) = CR(3,3);  
            
            STDCRA(j,1) = sqrt(CR(1,1));
            STDCRB(j,1) = sqrt(CR(2,2));
            STDCRR1(j,1) = sqrt(CR(3,3));
            
            T1(j,1) = 1/ML(3,1);
            T1LS(j,1) = 1/LS(3,1);
            STDCRT1(j,1) = 1*sqrt(CR(3,3))/(ML(3,1)^2);
            
            
            %dlx=nonzeros(EigVec(1,1,j)).*sqrt(nonzeros(EigVal(1,1,j)));
            %dux=nonzeros(EigVec(1,2,j)).*sqrt(nonzeros(EigVal(2,2,j)));
            %dly=nonzeros(EigVec(2,1,j)).*sqrt(nonzeros(EigVal(1,1,j)));
            %duy=nonzeros(EigVec(2,2,j)).*sqrt(nonzeros(EigVal(2,2,j)));
            
            %if(EigVal(1,1,j)>EigVal(2,2,j))
            %    deriv(j)=dly/dlx;
            %else
            %    deriv(j)=duy/dux;
            %end
            
 end  

%lx=nonzeros(EigVec(1,1,:)).*sqrt(nonzeros(EigVal(1,1,:)));
%ux=nonzeros(EigVec(1,2,:)).*sqrt(nonzeros(EigVal(2,2,:)));
%ly=nonzeros(EigVec(2,1,:)).*sqrt(nonzeros(EigVal(1,1,:)));
%uy=nonzeros(EigVec(2,2,:)).*sqrt(nonzeros(EigVal(2,2,:)));

snr=zeros(nrclusters,1);

for l=0:(nrclusters-1)   
    meanr2=mean(R1((l*nroccur)+1:(l+1)*nroccur));
    stdr2=std(R1((l*nroccur)+1:(l+1)*nroccur));
    snr(l+1,1)=meanr2/stdr2;
    %snr(l+1,2)=R1real(l*nroccur);
end

snrR1=R1realmore(1,:);
snrA=Bs';
snr=reshape(snr,nrR1s,nrBs);

%To show an A vs. R1 map with errorbars(either parallel with x and y axis
%or in the direction of the eigenvectors of the covariancematrix) enable one of the following lines: 

% errorbarxy(nonzeros(B)',nonzeros(R1)',sqrt(nonzeros(CRB))',(sqrt(nonzeros(CRR1)))',[],[],'w')
% %errorbarxyoblique(nonzeros(A)',nonzeros(R1)',lx',ly',ux',uy','w')
% hold on
% scatter(nonzeros(B),nonzeros(R1),'.k')


 scatter(reshape(repmat(TI,1,nrpoints),sizex*sizey,1),reshape(image,sizex*sizey,1),'.b')
 hold on;
 ti=0:1:round(max(TI));
%  yy=median(A)-median(B)*exp(-ti*median(R1)); 
%  yyls=median(Als)-median(Bls)*exp(-ti*median(R1ls));
%  yytrue=abs(As - Bs*exp(-ti*R1real(1)));
 yy=median(A)*(1-median(B)*exp(-ti*median(R1))+exp(-(TR(1)+ti)*median(R1)));
 yyls=median(Als)*(1-median(Bls)*exp(-ti*median(R1ls))+exp(-(TR(1)+ti)*median(R1ls)));
 yytrue=abs(As*(1-Bs*exp(-ti*R1real(1))+exp(-(TR(1)+ti)*R1real(1))));
 
%  plot(1:1000,1/0.003)
%  hold on
%  boxplot([T1, T1LS])
%  hold on
 
 plot(ti,abs(yy),'k');
 hold on
 plot(ti,yytrue,'r');
 hold on
 axis([0 max(ti)+100  min(min(image))-100 max(max(image))+250])
 plot(ti,yyls,'y')
%enable this line to create a surface plot of the snr a function of A and R
%surf (snrA,snrR1,snr, 'DisplayName', 'snrA'); 


