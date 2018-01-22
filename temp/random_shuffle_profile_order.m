% shuffle profile for random order 

for i=1:5
    R=randperm(size(profile_order,2));
profile_order_random(:,:,i)=profile_order(:,R,i);
end% shuffle 

profile_order=profile_order_random;
%%
filename='profile_22_1_2018_random_ordering'

savemask_LRT(profile_order,filename,visualize,modify_NSA_option)
fprintf('Saved as %s \n',filename)
fprintf('in folder: %s \n',pwd)
fprintf('Finished! \n \n \n ')
