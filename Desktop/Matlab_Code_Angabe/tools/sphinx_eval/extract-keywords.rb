#!/usr/bin/env ruby
# extract keywords for chime do_score script
# example: extract_keywords.rb result/transcription.txt

NUMBERS = {'ZERO' => 'z', 'ONE' => 1, 'TWO' => 2, 'THREE' => 3, 'FOUR' => 4, 'FIVE' => 5, 'SIX' => 6, 'SEVEN' => 7, 'EIGHT' => 8, 'NINE' => 9}
def word_to_number(word)
  NUMBERS[word.upcase] || '?'
end

File.read(ARGV[0]).each_line do |l|
  parts = l.split
  filename = parts[-2].match(/([a-z0-9_]+)$/)[1]
  if parts.size >= 7
    puts "#{filename} #{parts[3].downcase} #{word_to_number(parts[4])}"
  else
    # recognition failed; make a "random" guess
    puts "#{filename} s 2"
  end
end
