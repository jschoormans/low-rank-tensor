function MSE=visualize_convergence(iter,MSE,G,C,Phi,I,sdu,x,y)
% iter: iteration
% MSE: (first iter use [])
%I gold standard image
% sdu: size of tensor [x y te b]
% x,y": coordinates of pixel to track

    current_guess=G*C*Phi;
    MSE(iter)=sqrt(sum(abs(current_guess(:)-I(:)).^2));   
    
    cgr=abs(reshape(current_guess,sdu));
    
    figure(999);
    
    subplot(221);
    hold on 
    imshow(abs(cgr(:,:,1,1)),[]);
    plot(x,y,'r.','MarkerSize',20)
    hold off
    title('image at current iteration')
    
    subplot(223);
    hold on 
    imshow(abs(I(:,:,1,1)),[]);
    plot(x,y,'r.','MarkerSize',20)
    hold off
    title('Gold standard')
    
    subplot(322)
    plot([1:iter],MSE); title(MSE)
    xlabel('iter')
    title('Mean squared error')
    
    subplot(324)
    hold on 
    plot(squeeze(cgr(x,y,:,1)),'g')
    plot(squeeze(I(x,y,:,1)),'--k')
    hold off
    title(['DIM3: pixel value of x=', num2str(x),' y=',num2str(y)])
    
    subplot(326)
    hold on 
    plot(squeeze(cgr(x,y,:,1)),'g')
    plot(squeeze(I(x,y,:,1)),'--k')
    hold off
    title(['DIM4: pixel value of x=', num2str(x),' y=',num2str(y)])

drawnow;
end