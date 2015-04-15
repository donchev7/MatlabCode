clean = wavread('speech.wav');
noisy = clean + 0.01*randn(size(clean));
calcpesq(clean, noisy)