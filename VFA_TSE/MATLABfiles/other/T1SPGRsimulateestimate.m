%function [A,CRA,R2,CRR2] =simestimate(image, TE, sigma)

%simulate an A and R1 SPGR estimating experiment
%Henk Smit
%12.01.2011

clear('all')

%the 5 parameters you can change for a simulation, with in [] the value they
%have now

Angle=pi*(2:1:85)./180; %the echotimes in ms [5,15,..55]
Angle=Angle';
Es=[0.95:0.01:0.95]; %the R2's of the simdata [10,20...70]
As=[20:200:20];%the A's of the simdata [200,300..900]
Esmore=Es(ones(1,500),:); %the size of each A/R2 cluster(number of repetitions) [50]
sigma=15; %the sigma of the noise in the data [9]

Es=reshape(Esmore,1,numel(Esmore));

nroccur=size(Esmore,1);%
nrR2s=size(Es,2)/nroccur;
nrAngles=size(Angle,1);
nrAs=size(As,2);
nrpointsperR2=nroccur*nrR2s;
nrclusters=nrAs*nrR2s;
nrpoints=nrAs*nrR2s*nroccur;


Amat=As(ones(1,nrAngles),:);

for g=0:(size(As,2)-1)
    image(:,1+(g*nrpointsperR2):(g+1)*nrpointsperR2)=Amat(:,(g+1)*ones(1,nrpointsperR2)).*(sin(Angle)*(1-Es))./(1-(cos(Angle)*Es));
end

%Noise can be added by enabling the following three lines
noise=sigma*(randn(size(image,1),size(image,2))+i*randn(size(image,1),size(image,2)));
image=image+noise;
image=abs(image);

sizey=size(image,2);
sizex=size(image,1);

A=zeros(sizey,1);
CRA=A;
E=A;
CRE=A;
Els=A;
Als=A;

opt = optimset('fminunc');
opt = optimset(opt,'DerivativeCheck', 'off','Diagnostics','off','LargeScale','off','gradObj','on','Display','off','MaxIter',50,'Hessian','off','TolFun',1e-19,'Tolx',1e-9);
s = fitoptions('Method','NonLinearLeastSquares','Lower',[0,0,],'Upper',[30000,30000],'Startpoint',[500, 0.9],'Maxiter',20);
f = fittype('a*(1-b)*sin(x)/(1-b*cos(x))','options',s);

for j=1:sizey
    disp(['Calculating pixel: ' num2str(j) '/' num2str(sizey)]);
          for k=1:sizex
             yd(k)=double(image(k,j));
          end

            [LS,ML]=T1SPGRComputePar(yd',Angle,0,sigma,max(yd),opt,s,f);
            [CR]=T1SPGRCramerRao(ML,Angle,sigma);
            %[c2] = fit(echoes,data,f,s);
            %DT_LS(1,1)=c2.a;
            %DT_LS(2,1)=c2.b;
            %DT_LS(3,1)=c2.c;

            %[LS,ML,covmat]=ComputePar(yd',Angle);
            %[CR]=CramerRao(ML,Angle,sigma);
            %[EigVec(:,:,j),EigVal(:,:,j)]=eig(CR);
            
            E(j,1) = ML(2,1);
            A(j,1) = ML(1,1);
            Els(j,1)=LS(2,1);
            Als(j,1)=LS(1,1);
            CRA(j,1) = CR(1,1);
            CRE(j,1) = CR(2,2);  
            
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
    meanE=mean(E((l*nroccur)+1:(l+1)*nroccur));
    stdE=std(E((l*nroccur)+1:(l+1)*nroccur));
    snr(l+1,1)=meanE/stdE;
    %snr(l+1,2)=Es(l*nroccur);
end

snrR2=Esmore(1,:);
snrA=As';
snr=reshape(snr,nrR2s,nrAs);

%To show an A vs. R2 map with errorbars(either parallel with x and y axis
%or in the direction of the eigenvectors of the covariancematrix) enable one of the following lines: 

errorbarxy(nonzeros(A)',nonzeros(E)',sqrt(nonzeros(CRA))',(sqrt(nonzeros(CRE)))',[],[],'w')
%errorbarxyoblique(nonzeros(A)',nonzeros(R2)',lx',ly',ux',uy','w')
hold on
scatter(nonzeros(A),nonzeros(E),'.k')

%enable this line to create a surface plot of the snr a function of A and R
%surf (snrA,snrR2,snr, 'DisplayName', 'snrA'); 


