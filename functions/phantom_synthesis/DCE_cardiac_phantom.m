%
function I=DCE_cardiac_phantom(N,nframes,cardiac_period, varargin)
%simulates phantom ellipses (some get brighter), some pulsate with cardiac
%frequency.

%   N : resolution
%   frame: number of frames to simulate
%   cardiac period: (in frames) 
%   optional: visualize true/false

n_ellipse=6 %hard-coded for now...
%% EXAMPLE 
% to do...
%%

if nargin>3
    visualize = varargin{1};
    complexsim=varargin{2};
else
    complexsim=false;
    visualize =false;
end

if complexsim
    error('not yet implemented, sorry!'); end

%% START CODE

for nphantom=1:n_ellipse
    E(1)=0.4+0.6*rand; %intensity
    E2=rand*0.2; %length
    E3=rand*0.2; %width
    E(4)=-0.5+(rand); %x-coord of middle
    E(5)=-0.5+(rand); %y-coord of middle
    E(6)= rand*360; %angle
    for tt=1:nframes
        E(2)= 0.1+abs(E2*sin(pi*tt/cardiac_period)); %length
        E(3)= 0.1+abs(E3*sin(pi*tt/cardiac_period)); %width
        P{nphantom,tt} = phantom(E,N) ;
        
    end
end

% prepare DCE curves - very basic for now 
for nphantom=1:n_ellipse
    a=randn; 
DCE{nphantom}=a*(1+tanh(linspace(-3,rand()*10,nframes)));
end



% make simulation images
I=zeros(N,N,nframes);


if complexsim %% adds random phase...
  for nphantom=1:n_ellipse
    phase(nphantom)=exp(1i*rand*2*pi);
    fprintf('phase phantom %d is %4.2f \n',nphantom,angle(phase(nphantom)))
  end  
end


for ii=1:nframes
        for nphantom=1:n_ellipse % for all separate phantoms
%             if complexsim %add random phase to phantom
                % to do... 
%             else
                I(:,:,ii)=I(:,:,ii)+(P{nphantom,ii}).*DCE{nphantom}(ii);
%             end
        end
end

if visualize == 1
    % TO DO ...
    error('Not yet implemented, sorry!')
end

end

