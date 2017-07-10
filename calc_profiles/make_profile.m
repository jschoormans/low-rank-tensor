%calculate scan profile for LRT scan 
%where dimension 1 is the TSE echo number 

nDim1=60; % TSE dimensions
nDim2=5; % T2-prep 

ky=128; 
kz=128; 

ctrsize=5;
undersampling=0.01; %excluding centers

waiting_time=500e-3; 
TR=4e-3         %TR in ms; 
TR_shot=nDim1*TR+waiting_time; 
MC_maxiter=10000; 

nr_centerpoints=(2*ctrsize+1)^2; %number of k-points in the center squares; 
nr_points=ceil(undersampling*ky*kz); 
assert(nr_points>=nr_centerpoints,'fully sampled centers too big relative to undersampling')
nshots=nDim2*nr_points; 
total_time= TR_shot*nshots; %total time im seconds; 

fprintf('number of shots: %d, TSE number: %d, total time: %d seconds \n',nshots,nDim1,total_time)

%%
% add random point to every independent k-space, s.t. all have equal # of
% points
mask=zeros(ky,kz,nDim1,nDim2); 
fprintf('Starting Monte Carlo simulation \n')
for dim1=1:nDim1
    for dim2=1:nDim2
        fprintf('Dim 1: %d, Dim 2: %d |',dim1,dim2)
        m=zeros(ky,kz); MC_niter=0;
        if dim1==1 || dim2==1
            m=addCtr(m,ctrsize);
            while sum(m(:))~=nr_points
                m=rand(size(m))>(1-((nr_points-nr_centerpoints)/(ky*kz)));
                m=addCtr(m,ctrsize);

                MC_niter=MC_niter+1;
                if MC_niter>MC_maxiter
                    error('too many MC iterations - check settings')
                end
                
            end
        else
            
            while sum(m(:))~=nr_points
                m=rand(size(m))>(1-undersampling);
                MC_niter=MC_niter+1;
                if MC_niter>MC_maxiter
                    error('too many MC iterations - check settings')
                end
            end

        end
        fprintf(' Monte Carlo iterations: %d \n',MC_niter)
        mask(:,:,dim1,dim2)=m;
    end
end
clear m
%% 
figure(1); 
imshow(reshape(permute(mask,[1 3 2 4]),[ky*dim1,kz*dim2]))

%% Calculate profile ordering 

vec= @(x) x(:); 
% idea: for every shot: measure  points that are close in ky,kz after each
% other: for all shots sort points as such 
% inter-shot: permute such that inter shot ordering is random 

for dim2=1:nDim2
    m=mask(:,:,:,dim2);
    for dim1=1:nDim1
        [row,col]=find(mask(:,:,dim1,dim2))
        row=row-floor((1+ky)/2);
        col=col-floor((1+kz)/2);

        % sort into a certain number of radial spokes: low to high 
        [theta, rho]=cart2pol(row,col);
        
%         [~,s_index_angle]=sort(theta);
        [~,s_index_angle]=sort(row);

        % take blocks of N angle pieces (plus a residual) 
        rho_sorted=col(s_index_angle); 
        
        N=15;
        number_of_pieces=floor(nr_points/N);
        residual_N=mod(nr_points,N);
        
        index_r=[];
        for ii=1:number_of_pieces
        [~,index_piece]=sort(rho_sorted(1+(ii-1)*N:ii*N));
        index_piece_transform=s_index_angle((ii-1)*N+index_piece);
        index_r=[index_r;index_piece_transform];
        end
        [~,index_piece]=sort(rho_sorted(1+number_of_pieces*N:end));
        index_piece_transform=s_index_angle((number_of_pieces*N)+index_piece);
        index_r=[index_r;index_piece_transform];

%         figure(10); plot(row(index_r),col(index_r))

        kp(1,:,dim1)=row(index_r);
        kp(2,:,dim1)=col(index_r);
        
        % sort every block of alpha cols form low to high in rows
        
        
    end
    
    profile_order(1,:,dim2)=vec(permute(kp(1,:,:),[1 3 2]));
    profile_order(2,:,dim2)=vec(permute(kp(2,:,:),[1 3 2]));
    
end




%%  temp
clf
figure(2) 

for i=1:60:20000; 
    hold on 
    plot(profile_order(1,1:i),profile_order(2,1:i),'k.')
    plot(profile_order(1,i:i+60),profile_order(2,i:i+60),'r.')
    hold off
    xlim([-64 64])   
    ylim([-64 64])
pause(0.2)
    drawnow;
end 

%% 

p=profile_order;
delta_x=(p(1,1:end-1)-p(1,2:end));
delta_y= (p(2,1:end-1)-p(2,2:end));
figure(3)
plot(sqrt((delta_x).^2+(delta_y).^2),'.')




