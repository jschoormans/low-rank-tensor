% try  diffusion and T1 phantom 
% makes image for different b values and echo times
% then this image tensor is described as a low-rank tensor 

% phantom parameters
res=128;

ADCvals=[0.1 0.009 1 3 6];
T1vals=[20 50 100 200 800].*1e-3;

TE=[10:10:100].*0.001 % in seconds
bvals=3*[0:100:900].*1e-3 %units?

I=diffusion_T1_phantom(res,ADCvals,T1vals,TE,bvals,1);
%
figure(1); Q=[]
for ii=1:size(I,3)
    J=[];
    for jj=1:size(I,4);
        J=[J,abs(I(:,:,ii,jj))];
    end
    Q=[Q;J];
end
imshow(abs(Q),[0 1])
xlabel('echo times')
ylabel('b-values')
%% TUCKER MODEL 
aa=5:3:65
for iter=1:length(aa)
    a=aa(iter)
b=a
c=2;
d=2;

[F1,F2,F3,F4]=tucker(I,[a,b,c,d]);

disp('percentage of original information needed:')
frac(iter)=(a*b*c*d)./(res.^2*size(I,3)*size(I,4))
disp('frct of orignal signal explained:' ); 
pct(iter)=F3./100; 

end
figure(2)
plot(frac,pct);
xlabel('fraction of original signal size')
ylabel('explained variance by LRT model')

%% make simple Tucker model - plot all answers 


a=30
b=30
c=3
d=3
[F1,F2,F3,F4]=tucker(I,[a,b,c,d]);

size(F2)
figure(3)
subplot(2,2,1);imshow(abs(F1{1}))
title('decomposition 1')
subplot(2,2,2);imshow(abs(F1{2}))
title('decomposition 2')
subplot(2,2,3);imshow(abs(F1{3}))
title('decomposition 3')
subplot(2,2,4);imshow(abs(F1{4}))
title('decomposition 4')

figure(4);Q=[]
for ii=1:size(F2,3)
    J=[];
    for jj=1:size(F2,4);
        J=[J,abs(F2(:,:,ii,jj))];
    end
    Q=[Q;J];
end
imshow(abs(Q),[-1 1])
title('montage of core tensor'); colormap('jet')


% only image
irestore=zeros(res,res);
for ii=1
    for jj=1
   irestore= irestore+F1{1}*F2(:,:,ii,jj)*F1{2}.'     
    end
end

figure(5); 
imshow(abs(irestore),[])
title('F1{1}*C*F1{2}, first dimensions of core tensor are used')

 


