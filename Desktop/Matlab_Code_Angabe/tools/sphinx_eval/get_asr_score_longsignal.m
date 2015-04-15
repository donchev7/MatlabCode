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
function [ score ] = get_asr_score_longsignal(long_signal, listname, endpoints)
files = importdata([listname '.fileids']);
wavdir = tempname();
warning('off','MATLAB:MKDIR:DirectoryExists');
warning('off','MATLAB:audiovideo:wavwrite:functionToBeRemoved');
for file_ix=1:length(files)
    if file_ix==1
        startpoint=1;
    else
        startpoint=endpoints(file_ix-1)+1;
    end
    endpoint = endpoints(file_ix);
    mkdir(fileparts([wavdir '/' files{file_ix}]));
    %audiowrite([wavdir '/' files{file_ix} '.wav'], long_signal(startpoint:endpoint)./max(abs(long_signal)+0.01), 16000);
    wavwrite(long_signal(startpoint:endpoint)./max(abs(long_signal)+0.01), 16000, [wavdir '/' files{file_ix} '.wav']);
end
score=get_asr_score(wavdir, listname);
system(['rm -rf ' wavdir]);
