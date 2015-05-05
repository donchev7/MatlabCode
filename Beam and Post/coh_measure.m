function [varargout] = coh_measure(noise,ch1,ch2,alpha,mics,tit,N,L,fs)
% Measurement of the coherence function between two channels of a recorded multi-channel noise
% signal and comparison with the ideal sin(x)/x - coherence function
% [coh,coh_smooth] = coh_measure(noise,ch1,ch2,alpha,mics,tit,N,L,fs)
% coh
% coh_smooth deliveres the measured coherence function
% delivers a smoothed version of the meas. coherence
% function
% noise
% ch1
% ch2
% alpha input noise matrix recorded in a room
% first channel
% second channel
% factor for the exponentially weighted
% Welch periodogram; default = 0.8
% microphone position matrix; default mics_8xh.mat
% title of the plot; default tit = ''
% FFT length; default = 512
% decimation factor; default = 4
% sampling frequency; default = 16000 (16 kHz)
% mics
% tit
% N
% L
% fs
% Syntax:
% Plot measured coherence function:
% coh_measure(noise,ch1,ch2,alpha,mics,tit,N,L,fs)
% Save the measured coherence function and the smoothed c.f. in vectors:
% [coh,coh_smooth] = coh_measure(noise,ch1,ch2,alpha,mics,tit,N,L,fs)
% functions required: calc_cross.m

if nargin<9 fs = 16000; end
if nargin<8 L = 4; end
if nargin<7 N = 512; end
if nargin<6 tit = ''; end
if nargin<5 load mics_8xh.mat; end
if nargin<4 alpha = 0.8; end
if nargin<3
    help coh_measure
    return;
end
[K,Dim] = size(mics);
% Check microphone position matrix
if (Dim < 2) | (K < 1)
    error('bad matrix of microphine coordinates');
end
if Dim == 2
    rn = [mics zeros(K,1)];
else
    rn = mics;
end
% Define Distance Matrix of the Array
xc = rn(:,1);
xc = xc(:,ones(K,1));
dxc = xc - xc.';
yc = rn(:,2);
yc = yc(:,ones(K,1));
dyc = yc - yc.';
if Dim == 2
    dR = sqrt(dxc.^2 + dyc.^2);
else
    zc = rn(:,3);
    zc = zc(:,ones(K,1));
    dzc = zc - zc.';
    dR = sqrt(dxc.^2 + dyc.^2 + dzc.^2);
end
% set paramters
N2 = N/2 + 1;
n2 = 1:N2;
% Calc. spectral cross- and auto power density vectors
[cpsd,psd_x,psd_y] = calc_cross(noise,ch1,ch2,alpha,N,L,fs);
% Calc. coherence function
nom = sqrt(psd_x.*psd_y);
coh = cpsd./nom;
% Calc. microphone distance
d = dR(ch1,ch2);
if nargout == 0
    % Plot coherence functions
    coh_est = sinc((2*fs/N*d/340).*(n2-1))./(1 + 0.01);
    h = linspace(0,fs/2,N2);
    figure,plot(h,real(coh),'--b')
    hold on
    plot(h,coh_est,'r')
    hold off
    legend('Messung',['Theorie'],3)
    xlabel('Frequenz [Hz]')
    ylabel(['Real(\Gamma)'])
else
    % Save coherence functions in vectors
    varargout{1} = coh;
    varargout{2} = smooth(coh,30,'rloess');
end