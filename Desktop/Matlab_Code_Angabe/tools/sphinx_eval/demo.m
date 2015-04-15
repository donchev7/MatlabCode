% Andreas Schwarz (schwarz@lnt.de), 2013
%
% Example for simple ASR accuracy evaluation using Pocketsphinx.
%
addpath('sphinx_eval');
load('test_short');
% do some processing with long_signal
long_signal = fftfilt(hamming(50), double(long_signal)); % low-pass filter as example
% make sure output is (approximately) aligned with original signal
delay = 25;
long_signal = [long_signal(delay+1:end); zeros(delay,1)];
get_asr_score_longsignal(long_signal, 'test_short', endpoints) % keyword accuracy (as defined in CHiME challenge)