function [PESQ, FWSEGSNR, ASR] = evaluateScores(sig,y,cfg)

addpath(genpath('tools/fwsegsnr'));
addpath('tools/sphinx_eval');
addpath(genpath('tools/fwsegsnr'));

PESQ = calcpesq(sig.s(1:cfg.fs*100,1), y(1:cfg.fs*100));

param_fwsegsnr.frame = 0.025; %frame length in ms
param_fwsegsnr.shift = 0.01; %hop or shift in ms
param_fwsegsnr.window = @hanning; %window definition, needed for spectrogram
param_fwsegsnr.numband = 23;

FWSEGSNR=fwsegsnr(y(1:length(sig.s(:,1))), sig.s(:,1), cfg.fs, param_fwsegsnr);


test=load('test_100');
ASR=get_asr_score_longsignal(y, 'test_100', test.endpoints);


end