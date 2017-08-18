res1=size(kspace,1);
res2=size(kspace,2);
tensorsize=size(kspace);
sens=squeeze(sens)
%%
sens_normalized=bsxfun(@rdivide,sens,(eps+sum(abs(sens),3))); 
F=MCFopClass;
set_MCFop_Params(F,(squeeze(sens_normalized)),[res1,res2],[tensorsize(4),tensorsize(5)]);
fprintf('abs sum normalized sens \n')
fprintf('original kspace \n')
%  original kspace
k2=F*(F'*kspace); 
k3=F*(F'*k2); 
k4=F*(F'*k3); 

fprintf('%d %d %d %d \n',mean(abs(kspace(:))),...
mean(abs(k2(:))),...
mean(abs(k3(:))),...
mean(abs(k4(:)))...
)
fprintf('normalized kspace \n')

% normalize k-space
P0=F'*kspace;
[kspace_s,scaling]= scaleksp(kspace,P0); % scale kspace to ensure consistency over params;

k2=F*(F'*kspace_s); 
k3=F*(F'*k2); 
k4=F*(F'*k3); 

fprintf('%d %d %d %d \n',mean(abs(kspace(:))),...
mean(abs(k2(:))),...
mean(abs(k3(:))),...
mean(abs(k4(:)))...
)
%% original sens
sens_normalized=sens; 
fprintf('original sens \n')

F=MCFopClass;
set_MCFop_Params(F,(squeeze(sens_normalized)),[res1,res2],[tensorsize(4),tensorsize(5)]);
%  original kspace
k2=F*(F'*kspace); 
k3=F*(F'*k2); 
k4=F*(F'*k3); 
disp('original kspace')
fprintf('%d %d %d %d \n',mean(abs(kspace(:))),...
mean(abs(k2(:))),...
mean(abs(k3(:))),...
mean(abs(k4(:)))...
)
% normalize k-space
P0=F'*kspace;
[kspace_s,scaling]= scaleksp(kspace,P0); % scale kspace to ensure consistency over params;
k2=F*(F'*kspace_s); 
k3=F*(F'*k2); 
k4=F*(F'*k3); 
disp('normalize kspace')
fprintf('%d %d %d %d \n',mean(abs(kspace(:))),...
mean(abs(k2(:))),...
mean(abs(k3(:))),...
mean(abs(k4(:)))...
)

%% sum of sq sens
fprintf('sum of sq sens \n')

sens_normalized=bsxfun(@rdivide,sens,(eps+sqrt(sum(abs(sens).^2,3)))); 
F=MCFopClass;
set_MCFop_Params(F,(squeeze(sens_normalized)),[res1,res2],[tensorsize(4),tensorsize(5)]);
%  original kspace
k2=F*(F'*kspace); 
k3=F*(F'*k2); 
k4=F*(F'*k3); 

disp('original kspace')
fprintf('%d %d %d %d \n',mean(abs(kspace(:))),...
mean(abs(k2(:))),...
mean(abs(k3(:))),...
mean(abs(k4(:)))...
)
% normalize k-space
P0=F'*kspace;
[kspace_s,scaling]= scaleksp(kspace,P0); % scale kspace to ensure consistency over params;
k2=F*(F'*kspace_s); 
k3=F*(F'*k2); 
k4=F*(F'*k3); 
disp('normalize kspace')

fprintf('%d %d %d %d \n',mean(abs(kspace(:))),...
mean(abs(k2(:))),...
mean(abs(k3(:))),...
mean(abs(k4(:)))...
)



%% ?
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\t2shuffling-support-master\src'))
fprintf('sum of sq sens \n')
sensvec=reshape(sens,[size(sens,1)*size(sens,2),size(sens,3)]);

sens1_mag = reshape(vecnorm(reshape(sensvec, [], size(sens,3)).'), [size(sens,1),size(sens,2)]);
sens_normalized = bsxfun(@rdivide, sens, sens1_mag);
sens_normalized(isnan(sens_normalized)) = 0;
sens_normalized=reshape(sens_normalized,[size(sens)]);
% sens_normalized=sens;

F=MCFopClass;
set_MCFop_Params(F,(squeeze(sens_normalized)),[res1,res2],[tensorsize(4),tensorsize(5)]);

F_for = @(x) fft2c(x);
F_adj = @(y) ifft2c(y);

S_adj = @(as) (sum(bsxfun(@times, conj(sens_normalized), as), 3));
S_for = @(a) bsxfun(@times, sens_normalized, a);

A_adj = @(y) S_adj(F_adj(((y))));
A_for = @(a) ((F_for(S_for(a))));

kspace_s=kspace;
%
F=MCFopClass;
set_MCFop_Params(F,sens_normalized,[res1,res2],[tensorsize(4),tensorsize(5)]);

k2=F*(F'*kspace_s); 
k3=F*(F'*k2); 
k4=F*(F'*k3); 
fprintf('MCFLASS F %d %d %d %d \n',...
mean(abs(kspace_s(:))),...
mean(abs(k2(:))),...
mean(abs(k3(:))),...
mean(abs(k4(:))))


k2=F_for(F_adj(kspace_s)); 
k3=F_for(F_adj(k2));  
k4=F_for(F_adj(k3)); 
fprintf('FASJ FFOR FUNCTIONS %d %d %d %d \n',...
mean(abs(kspace_s(:))),...
mean(abs(k2(:))),...
mean(abs(k3(:))),...
mean(abs(k4(:))))


k2=A_for(A_adj(kspace_s)); 
k3=A_for(A_adj(k2));  
k4=A_for(A_adj(k3)); 
fprintf('AADJ FUNCTIONS %d %d %d %d \n',...
mean(abs(kspace_s(:))),...
mean(abs(k2(:))),...
mean(abs(k3(:))),...
mean(abs(k4(:))))
%%
F=MCFopClass;
set_MCFop_Params(F,sens_normalized,[res1,res2],[tensorsize(4),tensorsize(5)]);

tic
for i=1:5
    
    k2=F*(F'*kspace_s);
end
toc

%%
im1=F_adj(kspace);
% size(im1)
im2=S_for(S_adj(im1));
im3=S_for(S_adj(im1));
im4=S_for(S_adj(im1));
im5=S_for(S_adj(im1));

fprintf('%d %d %d %d %d\n',mean(abs(im1(:))),mean(abs(im2(:))),mean(abs(im3(:))),mean(abs(im4(:))),mean(abs(im5(:))))

