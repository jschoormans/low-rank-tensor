load S_temp_subspace1
S=S.';
L=5

[left_1,eigen_1,right_1]=svd(S,'econ');
eigenvals=diag(eigen_1);  %output for evaluation

[E1,E2]=eig(S*S');
eigenvals=abs(diag(gather(E2)));
nav_estimate=E1(:,1:L);
%%
figure(2)
subplot(321); 
plot(abs(diag(E2))./max(abs(E2(:))),'+'); hold on; plot(abs(eigenvals)./max(abs(eigenvals(:)))); hold off

subplot(323); 
plot(abs(E1(:,1:5)))
title('abs eigenvectors')
subplot(324); 
plot(abs(left_1(:,1:5)))
title('abs left singular values')

subplot(325); 
plot(angle(E1(:,1:5)))
title('angle eigenvectors')
subplot(326); 
plot(angle(left_1(:,1:5)))
title('angle left singular values')
%%
figure(3)
subplot(321); 
plot((abs(E1(:,1:5))-abs(left_1(:,1:5))))
title('difference in abs-vals')
subplot(322); 
plot((angle(E1(:,1:5))-angle(left_1(:,1:5))))
title('angle difference')




subplot(323); 
plot(real(E1(:,1:5)))
title('real eigenvectors')
subplot(324); 
plot(real(left_1(:,1:5)))
title('real left singular values')
subplot(325); 
plot(imag(E1(:,1:5)))
title('imag eigenvectors')
subplot(326); 
plot(imag(left_1(:,1:5)))
title('imag left singular values')

