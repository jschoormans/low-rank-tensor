spatial_images=reshape(G,[sdu(1) sdu(2) size(G,2)]);

I1=[];
I2=[];
for i=1:size(G,2);
I1=    cat(2,I1,abs(spatial_images(:,:,i)));
I2=    cat(2,I2,angle(spatial_images(:,:,i)));
end


%% plot spatial images for Ak 

spatial_images_Ak=reshape(Psi'*Ak,[sdu(1) sdu(2) size(G,2)]);

I3=[];
for i=1:size(Ak,2);
I3=    cat(2,I3,abs(spatial_images_Ak(:,:,i))./(eps+max(max(abs(spatial_images_Ak(:,:,i))))));
end

figure(9999);
subplot(311)
imshow(I1,[])
title('abs spatial images G_k')


subplot(312)
imshow(I3,[0 1])
title('abs spatial images A')
subplot(313)
imshow(I2,[-pi pi])
title('ph spatial images G_k')
