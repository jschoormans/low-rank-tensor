 function MSE=visualize_convergence(iter,MSE,G,C,Phi,I,sdu,x,y,kspace_1,F,Psi,Ak)
% iter: iteration
% MSE: (first iter use [])
%I gold standard image
% sdu: size of tensor [x y te b]
% x,y": coordinates of pixel to track

    x_idx = y;
    y_idx = x;

    current_guess=G*(C*Phi);
    calc_l_norms(F,Psi,current_guess,kspace_1,sdu,iter)
    
    run fig9999_spatial_images.m; 
    
    
    if ~isempty(I); 
%   MSE(iter)=sqrt(sum(abs(current_guess(:)-I(:)).^2))./numel(I);
%   %mean-squared error
    MSE(iter)=norm(abs(current_guess(:)-I(:)),'fro')^2/norm(abs(I(:)),'fro')^2; % relative error
    
    else
        MSE=[];
    end
    cgr=(reshape(current_guess,sdu));
    
    figure(999);
    
    subplot(221);
    imshow(abs(cgr(:,:,1,1,1)),[]);
    hold on 
    plot(x,y,'r+','MarkerSize',10)
    hold off
    title('image at current iteration')
    
    if ~isempty(I)
        subplot(223);
        imshow(abs(I(:,:,1,1)),[0 max(abs(cgr(:)))]);
        hold on
        plot(x,y,'r+','MarkerSize',10)
        hold off
        title('Gold standard')
    else
        subplot(223);
        imshow(angle(cgr(:,:,1,1,1)),[-pi pi]);
        hold on
        plot(x,y,'r+','MarkerSize',10)
        hold off
        title('image at current iteration (end,end)')
    end
    
    subplot(322)
        if ~isempty(I)
    plot([1:iter],MSE,'.-b'); 
    xlabel('iter')
    title('relative error')
        end
        
    h1=subplot(324);
    plot(squeeze(abs(cgr(x_idx,y_idx,1,:,1))),'g')
    if ~isempty(I)
        hold on
        plot(abs(squeeze(I(x_idx,y_idx,:,1))),'--k')
        hold off
    end
    title(['DIM3: pixel value of x=', num2str(x),' y=',num2str(y)])
    
    h2=subplot(326);
    plot(squeeze(abs(cgr(x_idx,y_idx,1,1,:))),'g')
    if ~isempty(I)
        hold on
        plot(abs(squeeze(I(x_idx,y_idx,1,:))),'--k')
        hold off
    end
    title(['DIM4: pixel value of x=', num2str(x),' y=',num2str(y)])

drawnow;
 end

 %===================================
 
 function     calc_l_norms(F,Psi,current_guess,kspace_1,sdu,iter)

     %==== Calculate l2-norm
    current_guess_k=F*current_guess;
    idx_nonempty=find(abs(kspace_1)>0);
    l2norm=sqrt(sum(abs(current_guess_k(idx_nonempty)-kspace_1(idx_nonempty)).^2));
    
    %==== Calculate l1-norm
    l1norm=sum(sum(abs(Psi*current_guess)));
    
    fprintf('l2-norm: %4.2f | l1-norm: %4.2f \n',l2norm, l1norm)
    
    %=== visualize in figure
    figure(3000);
    subplot(221); cla; title('l2 norm - abs'); 
    plot(abs(current_guess_k(idx_nonempty(1:sdu(1)))),'r'); hold on; plot(abs(kspace_1(idx_nonempty(1:sdu(1)))),'k+');
    hold off; legend('recon k-line','measured k-line');
    
    subplot(222); cla; title('l2 norm - phase'); 
    plot(angle(current_guess_k(idx_nonempty(1:sdu(1)))),'r'); hold on; plot(angle(kspace_1(idx_nonempty(1:sdu(1)))),'k+');
    hold off; legend('recon k-line','measured k-line');
    
    subplot(223); hold on; plot(iter,l2norm,'k*'); title('l2-norm ');
    subplot(224); hold on; plot(iter,sum(l1norm(:)),'k*'); title('l1-norm ')


 end
 