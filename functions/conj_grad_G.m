function Gk=conj_grad_G(G,C,Ak,Bk,Y,Z,alpha,Psi,du,Phi)
% Gamma_k = Omega(F G C_k Phi)
% Gamma*_k(d) should go from measured data to estimate of G....
% Omega : row selection operator??
% Psi: Wavelet operator

% in paper; conjugate gradient, but this is a linear least squares??? 


id=Psi'*(Ak+Y./alpha); %size (id) = res^2 x Lg (size of G) 


% find estimate of G_x based on data : adjoint omega with fixed C and Phi
%P=G C Phi
%G= P Phi.'

zf=fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(du,1),2),[],1),[],2),1),2); %zero filled recon
zf=reshape(zf,[size(zf,1)*size(zf,2),size(zf,3)*size(zf,4)]);
G_hat=zf*inv(Phi)*inv(C); % this is where the problem is...
disp('how to find G_hat???')

sum(G_hat(:))

nom=G_hat+(alpha/2).*id;


PsiPsi=Psi*(Psi')*ones(size(id)); %% really at a loss here....
LambdaLambda=ones(size(id)); % no idea about this one...S

denom=LambdaLambda+(alpha/2).*PsiPsi; % TO DO...
denom=denom.^(-1);

disp('to do: denominator ')
Gk=denom.*nom;
end