 function MSE=visualize_convergence(iter,MSE,G,C,Phi,I,sdu,x,y)
% iter: iteration
% MSE: (first iter use [])
%I gold standard image
% sdu: size of tensor [x y te b]
% x,y": coordinates of pixel to track

    x_idx = y;
    y_idx = x;

    current_guess=G*C*Phi;
    if ~isempty(I); 
    MSE(iter)=sqrt(sum(abs(current_guess(:)-I(:)).^2))./numel(I);   
    else
        MSE=[];
    end
    cgr=abs(reshape(current_guess,sdu));
    
    figure(999);
    
    subplot(221);
    imshow(abs(cgr(:,:,1,1,1)),[]);
    hold on 
    plot(x,y,'r+','MarkerSize',10)
    hold off
    title('image at current iteration')
    
    if ~isempty(I)
    subplot(223);
    imshow(abs(I(:,:,1,1)),[]);
    hold on 
    plot(x,y,'r+','MarkerSize',10)
    hold off
    title('Gold standard')
    end
    
    subplot(322)
        if ~isempty(I)
    plot([1:iter],MSE,'.-b'); title(MSE)
    xlabel('iter')
    title('Mean squared error')
        end
        
    h1=subplot(324);
    plot(squeeze(cgr(x_idx,y_idx,1,:,1)),'g')
    if ~isempty(I)
        hold on
        plot(abs(squeeze(I(x_idx,y_idx,:,1))),'--k')
        hold off
    end
    title(['DIM3: pixel value of x=', num2str(x),' y=',num2str(y)])
    
    h2=subplot(326);
    plot(squeeze(cgr(x_idx,y_idx,1,1,:)),'g')
    if ~isempty(I)
        hold on
        plot(abs(squeeze(I(x_idx,y_idx,1,:))),'--k')
        hold off
    end
    title(['DIM4: pixel value of x=', num2str(x),' y=',num2str(y)])

drawnow;
end