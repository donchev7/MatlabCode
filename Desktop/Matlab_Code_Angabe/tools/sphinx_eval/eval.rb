#!/usr/bin/env ruby

dataset = ENV['DATASET'] #|| "../test/BSE_andreas_2011-04-04_1551-nosilence"
modelbase = ENV['MODEL'] #|| '../sphinx_bse_map/models/BSE_andreas_2011-04-04_1551.cd_cont_601_mllr-map'
files = ENV['FILES'] || '../test.fileids'

datasetname = (dataset + '_' + modelbase.split(/\//).last).gsub(/[^a-z0-9\-_\.]/i,'_')
`mkdir -p results/#{datasetname}`
model = modelbase
`/CLUSTERHOMES/schwarz/clusterlocal/sphinx-20120222/bin/pocketsphinx_batch -debug 0 -hmm #{model} -jsgf chime.gram -fwdtree yes -fwdflat yes -bestpath no -dict #{model}/chime.dic -adcin yes -ctl #{files} -adchdr 44 -cepext .wav -cepdir #{dataset} -hyp results/#{datasetname}/transcription.txt 2>&1 | grep -v INFO`
#`~/local/pocketsphinx/bin/pocketsphinx_batch -debug 0 -hmm #{model} -jsgf chime.gram -fwdtree yes -fwdflat yes -bestpath no -dict #{model}/chime.dic -adcin yes -ctl #{files} -adchdr 44 -cepext .wav -cepdir #{dataset} -hyp results/#{datasetname}/transcription.txt 2>&1 | grep -v INFO`
`./extract-keywords.rb results/#{datasetname}/transcription.txt > results/#{datasetname}/result.txt`
score = `awk -f genscore.awk results/#{datasetname}/result.txt`
File.open("results/#{datasetname}/score.txt", 'w').write(score)
