function [ score ] = get_asr_score(wavdir, listname, force_regenerate )
currdir = fileparts(mfilename('fullpath'));
resultdir = [wavdir '/asr'];
score_filename = [resultdir '/chime_score.txt'];
if ~exist(score_filename, 'file') || (exist('force_regenerate','var') && force_regenerate)
    system(['rm -rf ' resultdir]);
    mkdir(resultdir);
    
    model = 'clean.cd_cont_600';
    
    % generate score
    [status, output] = system([currdir '/pocketsphinx_batch-0.8 -debug 0 -hmm ' fullfile(currdir, model) ' -jsgf ' currdir '/chime.gram -fwdtree yes -fwdflat yes -bestpath no -dict ' fullfile(currdir,model) '/chime.dic -adcin yes -ctl ' fullfile(currdir, listname) '.fileids -adchdr 44 -cepext .wav -cepdir ' wavdir ' -hyp ' resultdir '/transcription.txt']);
    if status ~= 0
        error('ASR failed. Output: %s', output);
    end
    system([currdir '/extract-keywords.rb ' resultdir '/transcription.txt > ' resultdir '/result.txt']);
    %[status, score_old] = system(['awk -f sphinx_eval/genscore.awk ' resultdir '/result.txt']);
    files = importdata([currdir '/' listname '.fileids']);
    result = fileread([resultdir '/result.txt']);
    if isempty(result)
       warning('ASR result file is empty %s:', wavdir); 
       score = NaN;
       return
    end
    lines = regexp(strtrim(result),'\s*\n\s*','split');
    if length(lines) < length(files)
        warning('only %d ASR results, expected %d', length(lines), length(files))
        score = NaN;
       return
    end
    score = calc_chime_score(lines);
    filewrite(score_filename, sprintf('%.2f\n',score));
    grid_score = calc_grid_score(fileread([resultdir '/transcription.txt']));
    filewrite([resultdir '/grid_score.txt'], sprintf('%.2f\n',grid_score));
else
    score = str2double(fileread(score_filename));
end
