function [d, e] = fwsegsnr(x, y, fs, param);
%% FWSNRSEG
%% Frequency-weighted segmental SNR
%%
%% [D, E] = FWSEGSNR(X, Y, FS, PARAM) calculates frequency-weighted segmental SNR of X
%% with reference to Y.
%%
%% Written and distributed by the REVERB challenge organizers on 1 July, 2013
%% Inquiries to the challenge organizers (REVERB-challenge@lab.ntt.co.jp)



% Normalization
%----------------------------------------------------------------------

x = x / sqrt(sum(x.^2));
y = y / sqrt(sum(y.^2));


% STFT
%----------------------------------------------------------------------

frame    = fix(param.frame * fs);
shift    = fix(param.shift * fs);
win      = window(param.window, frame);
noverlap = frame - shift;
fftpt    = 2^nextpow2(frame);

X = spectrogram(x, win, noverlap, fftpt, fs);
Y = spectrogram(y, win, noverlap, fftpt, fs);

X = abs(X);
Y = abs(Y);

[num_freq, num_frame] = size(X);


% Mel-scale frequency warping
%----------------------------------------------------------------------

melmat = fft2melmx(fftpt, fs, param.numband, 1, 0, fs / 2, 1, 1);
melmat = melmat(:, 1 : num_freq);

X = melmat * X;
Y = melmat * Y;


% Calculate SNR.
%----------------------------------------------------------------------

W = power(Y, 0.2);
E = X - Y;

ds = 10 * sum(W .* log10((Y.^2) ./ (E.^2)), 1) ./ sum(W, 1);
ds = min(ds, 35);
ds = max(ds, -10);


d = mean(ds);
e = median(ds);

