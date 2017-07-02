 function MSE=visualize_convergence(iter,MSE,G,C,Phi,I,sdu,x,y)
% iter: iteration
% MSE: (first iter use [])
%I gold standard image
% sdu: size of tensor [x y te b]
% x,y": coordinates of pixel to track

    x_idx = y;
    y_idx = x;

    current_guess=G*C*Phi;
    MSE(iter)=sqrt(sum(abs(current_guess(:)-I(:)).^2));   
    
    cgr=abs(reshape(current_guess,sdu));
    
    figure(999);
    
    subplot(221);
    imshow(abs(cgr(:,:,1,1)),[]);
    hold on 
    plot(x,y,'r+','MarkerSize',10)
    hold off
    title('image at current iteration')
    
    subplot(223);
    imshow(abs(I(:,:,1,1)),[]);
    hold on 
    plot(x,y,'r+','MarkerSize',10)
    hold off
    title('Gold standard')
    
    subplot(322)
    plot([1:iter],MSE,'.-b'); title(MSE)
    xlabel('iter')
    title('Mean squared error')
    
    h1=subplot(324);
    plot(squeeze(cgr(x_idx,y_idx,:,1)),'g')
    hold on 
    plot(squeeze(I(x_idx,y_idx,:,1)),'--k')
    hold off
    title(['DIM3: pixel value of x=', num2str(x),' y=',num2str(y)])
    
    h2=subplot(326);
    plot(squeeze(cgr(x_idx,y_idx,1,:)),'g')
    hold on 
    plot(squeeze(I(x_idx,y_idx,1,:)),'--k')
    hold off
    title(['DIM4: pixel value of x=', num2str(x),' y=',num2str(y)])

drawnow;
end