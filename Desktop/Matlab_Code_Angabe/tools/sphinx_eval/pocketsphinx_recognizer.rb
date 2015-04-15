require 'fileutils'

class PocketsphinxRecognizer
  def initialize(hmm=nil)
    # bad: "models/voxforge-en-0.3/model_parameters/voxforge_en_sphinx.cd_cont_3000"
    @hmm = hmm || "models/en_broadcastnews_16k_ptm256_5000"
  end

  def recognize(filenames, grammarfile, toprule=nil)
    File.open('/tmp/sphinx-file.txt', 'w+') do |f|
      filenames.each do |filename|
        f.puts(filename)
      end
    end
    output_dir = '/tmp/sphinx-out/' + rand.to_s
    FileUtils.mkdir_p(output_dir)
    model_dir = '.'
    puts "running speech recognizer..."
    system "/CLUSTERHOMES/schwarz/clusterlocal/sphinx-20120222/bin/pocketsphinx_batch -hmm #{@hmm} -jsgf #{grammarfile} #{toprule ? ("-toprule " + toprule) : ""} -fwdtree yes -fwdflat yes -bestpath yes -dict #{model_dir}/dicit.dic -adcin yes -ctl /tmp/sphinx-file.txt -adchdr 44 -cepext '' -cepdir '' -hyp #{output_dir}/transcription.txt"
    #system "~/local/pocketsphinx-0.7/bin/pocketsphinx_batch -hmm #{@hmm} -jsgf #{grammarfile} #{toprule ? ("-toprule " + toprule) : ""} -fwdtree yes -fwdflat yes -bestpath yes -dict #{model_dir}/dicit.dic -adcin yes -ctl /tmp/sphinx-file.txt -adchdr 44 -cepext '' -cepdir '' -hyp #{output_dir}/transcription.txt"
    puts "done."
    out = File.read(output_dir + '/transcription.txt')
    texts = out.lines.map do |l|
      l.match(/(.+?)\(/)[1].strip
    end
    if filenames.is_a?(String)
      return texts[0]
    else
      return texts
    end
  end
end
