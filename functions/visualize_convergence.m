 function MSE=visualize_convergence(iter,MSE,G,C,Phi,I,sdu,x,y)
% iter: iteration
% MSE: (first iter use [])
%I gold standard image
% sdu: size of tensor [x y te b]
% x,y": coordinates of pixel to track

%     x_idx = y; % ?
%     y_idx = x; % ?
    x_idx=x;
    y_idx=y;

    current_guess=G*C*Phi;
    run fig9999_spatial_images.m
    
    if ~isempty(I); 
%     MSE(iter)=sqrt(sum(abs(current_guess(:)-I(:)).^2))./numel(I);
%     %mean-squared error
    MSE(iter)=norm(abs(current_guess(:)-I(:)),'fro')^2/norm(abs(I(:)),'fro')^2; % relativer error
    
    else
        MSE=[];
    end
    cgr=(reshape(current_guess,sdu));
    
    figure(999);
    
    subplot(221);
    imshow(abs(cgr(:,:,1,8,2)),[]);
    hold on 
    plot(x,y,'r+','MarkerSize',10)
    hold off
    title('abs image (P) at current iteration')
    
    if ~isempty(I)
        subplot(223);
        imshow(abs(I(:,:,1,1)),[0 max(abs(cgr(:)))]);
        hold on
        plot(x,y,'r+','MarkerSize',10)
        hold off
        title('Gold standard')
    else
        subplot(223);
        imshow(angle(cgr(:,:,1,8,2)),[-pi pi]);
        hold on
        plot(x,y,'r+','MarkerSize',10)
        hold off
        title('angle image (P) at current iteration (end,end)')
    end
    
    subplot(322)
        if ~isempty(I)
    plot([1:iter],MSE,'.-b'); 
    xlabel('iter')
    title('relative error')
        end
        
    h1=subplot(324);
    plot(squeeze(angle(cgr(x_idx,y_idx,1,:,2))),'g')
    if ~isempty(I)
        hold on
        plot(angle(squeeze(I(x_idx,y_idx,:,2))),'--k')
        hold off
    end
    title(['DIM3: phase pixel value of x=', num2str(x),' y=',num2str(y)])
    
    h2=subplot(326);
    plot(squeeze(angle(cgr(x_idx,y_idx,1,7,:))),'g')
    if ~isempty(I)
        hold on
        plot(angle(squeeze(I(x_idx,y_idx,7,:))),'--k')
        hold off
    end
    title(['DIM4: phase pixel value of x=', num2str(x),' y=',num2str(y)])

drawnow;
end