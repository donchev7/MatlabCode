% unlike calc_chime_score, this function considers all words in the utterance
function [ relative_score ] = calc_grid_score( transcription )
correct = 0;
total = 0;

lines = regexp(strtrim(transcription),'\s*\n\s*','split');

for i=1:length(lines)
   line = lines{i};
   parts = regexp(line,'\(','split');
   recognized_words = regexp(strtrim(parts{1}),'\s+','split');
   filename = regexp(parts{2},'[A-Za-z0-9]{6}','match');
   filename=filename{1};
   for word_ix=1:length(filename)
      if word_ix > length(recognized_words) || isempty(recognized_words{word_ix})
          % wrong
      elseif lower(filename(word_ix)) == lower(recognized_words{word_ix}(1))
          correct = correct + 1;
      end
      total = total + 1;
   end
end

relative_score = 100*correct/total;

end

