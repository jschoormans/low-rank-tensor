function [M0, B0] = T2star_spins( T2star, nsamp , B0c, M0c)
% [M, B] = T2star_spins( T2star, nsamp [,B0, M0] );
% Creates a vector with spin magnitudes and offresonance 
% frequencies that is well suited for bloch simulations involving T2star.
%
% INPUTS:
%  T2star : T2star time of the simulated voxel(s), in ms
%  nsamp  : number of spins with which the T2star effect is simulated
%           Need many, default = 1001;
%  B0     : Base off-resonance frequency in rad/ms
%  M0     : Resting state Z-magnetisation of the voxel
% 
% OUTPUTS:
%  M      : fractional resting state Z-magnetisation for each spin
%  B      : Off resonance frequency for each spin.
%
% Use separate bloch simulation for each spin. 
% Then get the combined voxel value by summing the magnetisation vectors.
%
% Created by Dirk Poot, Erasmus MC,
% 25-1-2013

if nargin<2
    nsamp = 1001;
end;

xsmpls = (.5:nsamp)'/nsamp;
xsmpls_halfw = (1:nsamp-1)'/nsamp;
cxsmpls = cos(pi*xsmpls);
cxsmpls_halfw = cos(pi*xsmpls_halfw);
if nargin<4
    M0c = ones(size(T2star));
end;
cumpdf = [0; -1/pi*atan(tan( cxsmpls_halfw * pi/2))+.5; 1];
M0 = kron(M0c, diff(cumpdf)); 
B0 = kron(1./T2star , tan( cxsmpls * pi/2));
if nargin>=3 
    B0 = B0 + kron(B0c, ones(size(cxsmpls))) ; 
end;