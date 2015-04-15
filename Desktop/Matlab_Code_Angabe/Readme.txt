Short introduction to this code
* Main file: main.m
  Sets parameters with function SetParams() and creates microphone signals with LoadMicInputs()

* Folgende Parameter sind jetzt erstmal wichtig für dich (ist auch alles in SetParams.m kommentiert
  - cfg.nmic 		->	Number of microphones
  - cfg.RIRcond		->	Determines used set of HRIRs
  - cfg.path_source 	->	Path to source signals
  - cfg.mic_ch		->	Microphone channels. Length(cfg.mic_ch) has to be consistent with cfg.nmic				
  - cfg.sig_len		->	Length of signals in seconds
  - cfg.noise_type	->	In order to add additional (self- or background) noise to the microphone signals
