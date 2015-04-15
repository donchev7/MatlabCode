function snr = MeasureAliasDist(h,K,N);
% snr = MeasureAliasDist(h,K,N);
%
% For white gaussian noise input, calculate an SNR measure for the resulting
% aliasing in an oversampled GDFT filterbank with K channels, subsampled by N.
%
% Input parameters:
%      K        twice number of subband channels used
%      N        decimation ratio in the subbands
%
% Output parameter:
%      snr      signal-to-noise ratio due to aliasing.
%
% St.Weiss, University of Strathclyde, 18.8.1997

Lh = length(h);
% * old *n = Lh*N;           % number of frequency points measured.
n = Lh;
H = freqz(h,1,n);

% *old* Baseband = sum(abs(H(1:Lh)).^2); 
% *old* Aliasing =  sum(abs(H(Lh+1:Lh*N)).^2);
Baseband = sum(abs(H(1:floor(n/N))).^2); 
Aliasing =  sum(abs(H(ceil(n/N)+1:n)).^2);

snr = 10*log10(Baseband/Aliasing);

