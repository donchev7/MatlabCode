function p = calcpesq(ref, target)
%% CALCPESQ
%% PESQ based on ITU-T Recommendation P.862.2

oldpath = pwd();
tmppath = tempname();
mkdir(tmppath);
cd(tmppath);

% Run PESQ software. 
%----------------------------------------------------------------------

if exist('pesq_results.txt', 'var')
  delete('pesq_results.txt');
end

scaling_factor = 0.999/max([max(abs(target)), max(abs(ref))]);


target = target(1:min([length(ref),length(target)]));
ref = ref(1:min([length(ref),length(target)]));

wavwrite(ref*scaling_factor,16000,'ref.wav');
wavwrite(target*scaling_factor,16000,'target.wav');

CMD = [fileparts(mfilename('fullpath')), '/pesq/PESQ +16000 +wb ref.wav target.wav > TMP'];
system(CMD);

% Extract the result.
%----------------------------------------------------------------------

fid = fopen('pesq_results.txt');

l = fgetl(fid);  %% first line is discarded
l = fgetl(fid);  %% second line contains the result

[a, l] = strtok(strtrim(l));  %% remove reference file name
[a, l] = strtok(strtrim(l));  %% remove target file name
[a, l] = strtok(strtrim(l));  %% remove target file name

p = str2num(strtok(strtrim(l)));

fclose(fid);

delete('pesq_results.txt');
delete('TMP');

cd(oldpath);
system(['rm -rf ' tmppath]);