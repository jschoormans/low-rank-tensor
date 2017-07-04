Mxy=Mxy_muscle;
%%
figure;  set(0,'DefaultAxesFontSize', 16);
i=sqrt(-1);
K_space=zeros(max(abs(Ky_order))*2+1,max(abs(Kz_order))*2+1);
for(nn=1:length(Ky_order))
    intensity=Mxy(ETL_nr(nn));
    if(AQ_enable(nn)==1)
       K_space(max(abs(Ky_order))+1+Ky_order(nn),max(abs(Kz_order))+1+Kz_order(nn))=intensity;
    end
end
 imagesc(K_space,[0 1]); colormap jet; colorbar;
 axis equal;
PSF=abs(fftshift(fft2(K_space)));

%% retime plot
figure('Position',[100 100, 800, 800])
K_space=zeros(max(abs(Ky_order))*2+1,max(abs(Kz_order))*2+1);
for(nn=1:length(Ky_order))
    intensity=Mxy(ETL_nr(nn));
    if(AQ_enable(nn)==1)
        K_space(max(abs(Ky_order))+1+Ky_order(nn),max(abs(Kz_order))+1+Kz_order(nn))=intensity;
        imagesc(K_space,[0 1]); colormap jet; colorbar;
        axis equal;
        drawnow;
        pause(0.01);
    end
end
    

  