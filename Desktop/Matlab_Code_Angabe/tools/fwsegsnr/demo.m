% Params from REVERB challenge
param_fwsegsnr = struct('frame'  , 0.025, ...
    'shift'  , 0.01, ...
    'window' , @hanning, ...
    'numband', 23);

clean = wavread('speech.wav');
noisy = clean + 0.01*randn(size(clean));
fwsegsnr(noisy, clean, 16000, param_fwsegsnr)