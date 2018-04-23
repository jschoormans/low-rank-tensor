function [Ak,lambda]=soft_thresh_A(G,Y,alpha,lambda_o,Psi,operatorsize,params)
% applies elementwise soft-thresholding operator
% as described in eqn B1

lambda=lambda_o; 

if params.iter>1 || params.autolambda==0
    T=(Psi*G)-Y./alpha;
    % l2 norm of T
    l2T=sqrt(sum(abs(T).^2,2));
    thr=(l2T>lambda/alpha);
else        %first iter and automu -
    fprintf('Automatically selecting lambda...\n')
    
    nthr=0;
    iter_autolambda=0;
    while (nthr<0.3)||(nthr>0.9)
        iter_autolambda=iter_autolambda+1;
        
        if nthr<0.3;
            lambda=lambda.*0.9;
        elseif nthr>0.9
            lambda=lambda.*1.1;
        end
        if iter_autolambda>100
            lambda =lambda_o;
            fprintf('automu failed \n')
        end
        
        T=(Psi*G)-Y./alpha;
        l2T=sqrt(sum(abs(T).^2,2));
        thr=(l2T>lambda/alpha);
        
        nthr=sum(thr(:))./numel(thr);
        if iter_autolambda>100
            break;
        end
    end
        fprintf('lambda= %d \n',lambda)

end







 

fprintf('soft-thresholding A: ')
fprintf('lambda/alpha: %1.2e  |',(lambda/alpha))
fprintf('thresholding %1.2f  percent \n',(100*(sum(~thr(:))./numel(thr))))

st=thr.*((l2T-(lambda/alpha))./(l2T+eps)); %soft thresholding
Ak=repmat(st,[1 size(T,2)]).*T;

%visualise Ak this iteration (mag and phase, first two).
fprintf('dimension 64 82 and rank 6 assumed here! A visualisation removed for now \n')
% Akres=reshape(Ak,[64 82 6]);
% figure(12321);subplot(241);imshow(squeeze(abs(Akres(:,:,1))),[]);colorbar;title('First Ak (mag)');
% subplot(242);imshow(squeeze(angle(Akres(:,:,1))),[]);colorbar;title('First Ak (phase)');
% subplot(243);imshow(squeeze(abs(Akres(:,:,2))),[]);colorbar;title('Second Ak (mag)');
% subplot(244);imshow(squeeze(angle(Akres(:,:,2))),[]);colorbar;title('Second Ak (phase)');

if params.visualization == 1;
    title_text = sprintf('s.t. A: %d data points thresholded (%1.2f%%).',sum(~thr(:)), (100*(sum(~thr(:))./numel(thr))));
    figure(998);subplot(221); spy(reshape(thr,[operatorsize]));title(title_text); % removed for now 15-11-2017, unindented 30-11-2017
end
end
