function [ relative_score ] = calc_chime_score( lines )
cnt = 0;
score = 0;

for i=1:length(lines)
   line = lines{i};
   if line(4) == line(8)
       score = score+1;
   end
   if line(5) == line(10)
       score = score+1;
   end
   cnt = cnt+2;
end

relative_score = 100*score/cnt;

end

