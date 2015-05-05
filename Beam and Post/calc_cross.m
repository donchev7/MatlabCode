function [CXX,CX1,CX2] = calc_cross(signal,ch1,ch2,alpha,N,L,fs)
% [CXX,CX1,CX2] = calc_cross(signal,ch1,ch2,alpha,N,L,fs)
% Calculate spectral cross- and auto power density vectors for 2 channels of a
% recorded multi-channel signal
% CXX
% CX1
% CX2 cross power density vector
% auto power density vector of channel 1
% cross power density vector of channel 2
% signal
% ch1
% ch2
% alpha input signal matrix recorded in a room
% channel 1
% channel 2
% factor for the exponentially weighted
% Welch periodogram; default = 0.8
% FFT length; default = 512
% decimation factor; default = 4
% sampling frequency; default = 16000 (16 kHz)
% N
% L
% fs
% used function: welch_est.m
if nargin<7 fs = 16000; end
if nargin<6 L = 4; end
if nargin<5 N = 512; end
if nargin<4 alpha = 0.8; end
if nargin<3
    help calc_cross
    return;
end
M = N/L;
% Zero-padding to reach a signallength to be a multiple of L
[Nx,K] = size(signal);
dum = ceil(Nx/N)*N - Nx;
signal = [signal;zeros(dum,K)];
Nx = length(signal);
% Initialise Vectors and Matrices
h = hanning(N);
H = h(:) * ones(1,K);
N2 = N/2 + 1;
y = zeros(Nx,K);
n2 = 1:(N2);
% initialise power density vectors
CXX = zeros(N2,1);
CX1 = zeros(N2,1);
CX2 = zeros(N2,1);
for k = 1:M:(Nx - 2*N + 1)
    k1 = k:k+N-1;
    % FFT - Filterbank with Hanning-Windowing
    X = fft(signal(k1,:) .* H,N).';
    % Calc. spectral power density vectors
    Y = (X(:,n2)).';
    CXX = welch_est(CXX,Y(:,ch1),Y(:,ch2),alpha);
    CX1 = welch_est(CX1,Y(:,ch1),Y(:,ch1),alpha);
    CX2 = welch_est(CX2,Y(:,ch2),Y(:,ch2),alpha);
end