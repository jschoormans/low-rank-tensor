guess=(F*(G*C*Phi.'));
du_shift=du;%ifftshift(ifftshift(du,1),2);
% du_shift=du;%ifftshift(ifftshift(du,1),2);

dur=reshape(du_shift,[128^2,100]);

guess(:)'*guess(:)
dur(:)'*dur(:)
imag(guess(:)-dur(:))'*imag(guess(:)-dur(:))
real(guess(:)-dur(:))'*real(guess(:)-dur(:))


im=F'*dur; 
im2=F'*(F*(F'*dur)); 

guessim=F'*guess; 
guessim2=F'*(F*(F'*guess)); 

close all
figure(1); hold on; plot(abs(guess(:,10))); plot(abs(dur(:,10))); hold off

figure(2);subplot(221);
imshow(reshape(abs(im(:,1)),[128 128]),[]);
subplot(222);
imshow(reshape(abs(im2(:,1)),[128 128]),[]);
subplot(223);
imshow(reshape(abs(guessim(:,1)),[128 128]),[]); 
subplot(224);
imshow(reshape(abs(guessim2(:,1)),[128 128]),[]); 


figure(3);subplot(221);
imshow(reshape(imag(im(:,12)),[128 128]),[]);
subplot(222);
imshow(reshape(abs(im2(:,12)),[128 128]),[]);
subplot(223);
imshow(reshape(imag(guessim(:,12)),[128 128]),[]); 
subplot(224);
imshow(reshape(abs(guessim2(:,12)),[128 128]),[]); 

figure(4); 
subplot(224);
imshow(reshape(abs(P0(:,1)),[128 128]),[]); % same as guess images... --> Fdu is wrong
