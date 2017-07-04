
        A=4000;
        B=1.6;
        T1=1200;
        R1=1/T1;
        T2=60;
        R2=1/T2;
        TI=[50 100 200 400 800 1600 2400 4000]';
        TR=[2000 2000 2000 2000 2000 2000 2000 2000]';
                initval = [0.85, 4050, 1.65]';


        E11 = exp(-TI*R1);
        E1 = bsxfun(@times,B,E11);
        E2 = exp(-TR*R1);
        S = bsxfun(@times, A, bsxfun(@plus, 1, bsxfun(@plus, E2 , -E1)));
        

for k=1:20
    for j=1:100



        S2=AddRiceNoise(S,10);

    %     plot(TI,S)
    %     hold on
    %     plot(TI,S2,'k')
        yd=S2;



        TRe=((k/20)+0.5)*[2000 2000 2000 2000 2000 2000 2000 2000]';

         fun = @(tht) predict_IRT1TR( tht , TI,TRe);
         tht0 = fit_MRI( fun, yd, initval); 
         tht4 = fit_MRI( fun, yd, tht0); 
         tht2 = fit_MRI( fun, yd, tht4); 
         tht3 = fit_MRI( fun, yd, tht2); 
         tht1 = fit_MRI( fun, yd, tht3); 

         T1es(j,k)=1000/tht1(1);

    %      disp(['T1= ' num2str(1000./tht1(1)) ' A= ' num2str(tht1(2)) ' B= ' num2str(tht1(3))])

    end
end

mean(T1es)

legend('True TR = 2000')













