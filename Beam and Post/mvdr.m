function [w,d0] = mvdr(tau,Gamma,f,nr)
% [w,d0] = mvdr(tau,Gamma,f,nr)
% Compute weights of MVDR-Beamformer
% w
% d0 Weight vector of length K (K ... number of mics)
% Steering vector
% tau
% Gamma
% f
% nr Time alignment vector
% Coherencematrix
% Frequency in Hz
% you can choose between the following beamformers
% 'DSB' ... Delay&Sum-Beamformer
% 'SDB' ... Superdirective Beamformer

if nargin <4
    help mvdr;
    return;
end
[K,dum] = size(Gamma);
beta = 2*pi*f/340;
% wave number
% calc. steering vector of desired direction
d0 = exp(-1i * beta * tau);
switch nr
    case 'DSB'
    % Delay-Sum-Beamformer
        w = d0 / K;
    case 'SDB'
        % Superdirective Beamformer for Diffuse Noise Field
        B = (Gamma^-1)*d0;
        Lambda = (d0'*B)\1;
        % Lagrange multiplicator
        w = B*Lambda;
        % optimum coefficient vector at given frequency
    otherwise
        error('Please insert correct number for the beamformer you want to choose!')
end