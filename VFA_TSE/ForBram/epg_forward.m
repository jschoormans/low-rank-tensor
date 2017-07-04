function [s,m_long,FZall] = epg_forward(flipangle,numberechos,esp_vec,T1,T2,m0_frac)
% Forward EPG algorithm for computing signal values
% Output:
%   s = signal (real-valued)
%   FZall = array containing F and Z states
%
% Input:
%   flipangle: FA array in degree  e.g. [180 120 58 37 25 10...]
%   numberechos: echo train length  e.g. 6
%   esp_vec: echo spacing vector in ms   e.g. [3 3 3 3 3 3]
%   T1: T1 value in ms 
%   T2: T2 value in ms
%   m0_frac: m0 at the begining of TSE sequence. Normally 1
% 
% 

if (isempty(numberechos))
    numberechos = length(flipangle);
elseif (length(flipangle)==1 && numberechos>1 && abs(flipangle)<pi)
    flipangle(2) = flipangle(1);
    flipangle(1) = (pi*exp(1i*angle(flipangle(2)))+flipangle(2))/2;
% (Henning's 1st flip reduced trick)
elseif (numberechos>length(flipangle))
    flipangle(end+1:numberechos) = flipangle(end);
end
% Convert degrees to radians:
flipangle = flipangle*pi/180;
% Allocate all known states:
FZ = zeros(3,2*numberechos);
% Prepare arrays to store all F and Z states:
FZall = cell(1,1,numberechos); % Add dephased states with each gradient application
% Prepare array to store signal:
s = zeros(1,numberechos);
m_long = zeros(1,numberechos); % longitudinal magnetisation
% Initial conditions after excitation:
FZ(1,1) = 1*m0_frac;
FZ(2,1) = 1*m0_frac;
% Refocusing pulses and dephasing:
for ech=1:numberechos
    FZ = epg_dephase(FZ,esp_vec(ech)/2,1,T1,T2); % Left crusher
    FZ = epg_rf(FZ,flipangle(ech)); % Refocusing RF pulse
    FZ = epg_dephase(FZ,esp_vec(ech)/2,1,T1,T2); % Right crusher
    % Store signal corresponding to current echo:
    s(ech) = abs(FZ(1,1));
    m_long(ech) = abs(FZ(3,1));
    % Store states:
    FZall{1,1,ech} = FZ;
end
end

function FpFnZ = epg_dephase(FpFnZ,t,Gon,T1,T2)
% Propagate EPG states through a period of dephasing due to relaxation and
% gradients. Describes time evolution between RF pulses
% Output:
%   FpFnZ = updated F+, F- and Z states

% Decay of states due to relaxation alone:
E1 = exp(-t/T1);
E2 = exp(-t/T2);
EE = diag([E2 E2 E1]);
RR = 1-E1; % Mz recovery, affects only Z0 state
% Relaxation (and recovery):
FpFnZ = EE * FpFnZ;
FpFnZ(3,1) = FpFnZ(3,1) + RR;
% Advance states:
if Gon
    FpFnZ = epg_grad(FpFnZ);
end
end

function FpFnZ = epg_rf(FpFnZ,alpha)
% Propagate EPG states through an RF pulse
% Input:
%   FpFnZ = 3xN vector of F+, F- and Z states
%   alpha = flip angle of pulse (radians)
% Output:
%   FpFnZ = Updated FpFnZ state

if (real(alpha) > 2*pi)
    warning('epg_rf:  Flip angle should be in radians!')
end
T = [(cos(alpha/2))^2 (sin(alpha/2))^2 sin(alpha);
    (sin(alpha/2))^2 (cos(alpha/2))^2 -sin(alpha);
    -0.5*sin(alpha) 0.5*sin(alpha) cos(alpha)];
% T mixes states of equal dephasing order:
FpFnZ = T * FpFnZ;
end

function FpFnZ = epg_grad(FpFnZ)
% Propagate states through a "unit" gradient
% Output:
%   FpFnZ = Updated FpFnZ state

% Gradient shifts dephasing order of F states:
FpFnZ = [FpFnZ zeros(3,20)]; % Add higher dephased states
FpFnZ(1,:) = circshift(FpFnZ(1,:),[0 1]); % Shift Fp states
FpFnZ(2,:) = circshift(FpFnZ(2,:),[0 -1]); % Shift Fn states
FpFnZ(2,end) = 0; % Zero highest Fn state
FpFnZ(1,1) = conj(FpFnZ(2,1)); % Fill in lowest Fp state
end