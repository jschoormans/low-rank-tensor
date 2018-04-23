function savemask_LRT(profile_order,varargin)
modify_NSA_option=0; %standard: 0: normal way of calculating NSA 
if nargin==1
    visualizeflag=1;
    filename='tempmask'
elseif nargin==2
    filename=varargin{1};
    visualizeflag=1;
elseif nargin==3
    filename=varargin{1};
    visualizeflag=varargin{2};
elseif nargin==4
    filename=varargin{1};
    visualizeflag=varargin{2};
    modify_NSA_option=varargin{3}; % if 1: recalculate NSA to prevent memory problems during scans
    
end


if ~strcmp(filename(end-3:end),'.dat') 
    filename=[filename,'.dat'];
end
    
profile_reshaped = reshape(profile_order, 2, size(profile_order, 2) * size(profile_order, 3))';


%%% append NSA to profiles 
nsa = zeros(size(profile_reshaped, 1), 1);

for k = 1 : size(profile_order, 3);
        profile_this_dynamic = profile_order(:,:,k)';
    for p = 1 : size(profile_order, 2);
        ky_match = (profile_this_dynamic(1:p-1, 1) == profile_this_dynamic(p, 1));
        kz_match = (profile_this_dynamic(1:p-1, 2) == profile_this_dynamic(p, 2));
        
        
        match_profiles_nr = sum( ky_match .* kz_match );
        nsa(p + (k-1) * size(profile_order, 2)) = match_profiles_nr;

        
        
    end
    
end
if modify_NSA_option; %ONLY FOR OFFLINE RECON- ONLINE RECON ON SCANNER SHOULD BE DISABLED!!!
    %calculate maximum number of NSA needed 
%     n_klines_per_dyn=size(profile_order,2) %number of lines per dyanmic  
%     n_lines_fully_sampled=(2*(max(profile_order(1,:)))+1)*(2*(max(profile_order(2,:)))+1); % ky*kz   
%     max_NSA=ceil(n_lines_fully_sampled/n_klines_per_dyn)-1 % -1 because you start counting at 0 ; 
    max_NSA=0; 
    nsa(nsa>max_NSA)=0; %delete the NSA's which are too much for memory allocation 
end


    profile_final=cat(2,profile_reshaped,nsa);
disp(['max nsa:',num2str(max(nsa))]);
if visualizeflag
    disp('to do: CHANGE CODE FOR VISUALIZATION')
% kspa = zeros(64, 128);
% for m = 1:9840;
%     x_idx = profile_this_dynamic(m, 1) + 64;
%     y_idx = profile_this_dynamic(m, 2) + 32;
%     kspa(x_idx, y_idx) = kspa(x_idx, y_idx) + 1;
% end
% figure(10); imshow(kspa,[]); colormap jet
end


dlmwrite(filename,profile_final,'delimiter','\t');
fprintf('profile saved as: %s \n',filename)
