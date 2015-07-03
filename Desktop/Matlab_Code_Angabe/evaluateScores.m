function [PESQ, FWSEGSNR, ASR] = evaluateScores(ref,y,cfg,sig)

addpath(genpath('tools/calcpesq'));
addpath('tools/sphinx_eval');
addpath(genpath('tools/fwsegsnr'));

PESQ = calcpesq(sig.s(1:cfg.fs*100,1), y(1:cfg.fs*100));
%PESQ = calcpesq(ref, y);

param_fwsegsnr.frame = 0.025; %frame length in ms
param_fwsegsnr.shift = 0.01; %hop or shift in ms
param_fwsegsnr.window = @hanning; %window definition, needed for spectrogram
param_fwsegsnr.numband = 23;

FWSEGSNR=fwsegsnr(y(1:length(ref(:,1))), ref(:,1), cfg.fs, param_fwsegsnr);
%FWSEGSNR=fwsegsnr(y, ref(1:length(y)), cfg.fs, param_fwsegsnr);

test=load('test_100');
test.listname = 'test_100';
ASR=get_asr_score_longsignal(y, 'test_100', test.endpoints);


end