function [lds] = welch_est(lds_prev,Xi,Xj,alpha)
% [korr] = welch_est(korr_prev,Xi,Xj,fs,D,tau)
% Calculates the Welch-Estimate of Auto- or Cross-Spectral-Density
% lds Auto- or Cross-Spectral Density for the current frame
% lds_prev
% Xi
% Yi Auto- or Cross-Spectral Density for the previous frame
% signal vector of signal i
% signal vector of signal j;
% default Xj = Xi (for Auto-Spectral density)
% Welch-factor for weighting the previous frame; default alpha = 0.8
% alpha
if nargin<4 alpha = 0.8; end
if nargin<3 Xj = Xi; end
if nargin<2
    help welch_est;
    return;
end
% Calc. recursive Welch-formula
first = alpha.*lds_prev;
delta = Xi.*conj(Xj);
delta = (1 - alpha).*delta;
% calc. the spectral density
lds = first + delta;