% get_asr_score_longsignal(long_signal, listname, endpoints)
%
% Andreas Schwarz (schwarz@lnt.de), 2013
%
% This function takes a concatenation of GRID utterances, splits it, calls the
% speech recognizer, and returns the keyword score (recognition rate for
% letter and number in the utterance, as in the CHIME challenge).
%
% Example:
%    addpath('sphinx_eval');
%    load('test_short');
%    % do something processing with long_signal; make sure output has same time alignment
%    long_signal = fftfilt(hamming(20), double(long_signal)); % low-pass filter as example
%    get_asr_score_longsignal(long_signal, 'test_short', endpoints)