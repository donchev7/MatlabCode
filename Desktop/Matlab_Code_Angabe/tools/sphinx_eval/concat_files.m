% This script takes a list of files, and concatenates the wave data of all
% files to a long vector. The end times of the utterances are saved in a
% variable.

% short: 100 utterances; medium: 500 utterances
%start index
start = 101;
%number of utterances to concatenate
number = 100;
%load list that includes all file ids
files=importdata(['test.fileids']); 
%set listname of list that will include the file ids
% of the utterances concatenated in the resulting test file
listname = ['test_short_' num2str(start) '_' num2str(start+number-1) '.fileids'];

%initialization
long_signal = [];
endpoints = [];
fileID = fopen(listname,'w');
if (fileID == -1) %error
   error('List of file ids could not be created!\n'); 
end
%concatenate corresponding utterances and create list of file ids
for file_ix=start:(start+number-1)
    signal=wavread(['/CLUSTERHOMES/LMS_AUDIO/schwarz/data-bse-dereverb/clean/' files{file_ix} '.wav']);
    long_signal = single([long_signal; signal]);
    endpoints(file_ix-number) = length(long_signal);
    %create list of file ids
    fprintf(fileID,[files{file_ix},'\n']);
end
%save files (you may need to rename these files according to your
%organization of data)
wavwrite(long_signal, 16000, ['test_' int2str(number) '.wav'])
save(['test_' int2str(number) '.mat'], 'endpoints', 'listname', 'long_signal');
%close list
fclose(fileID);