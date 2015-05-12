function beampattern(steerV,Wmvdr,Wrfsb,angRangeAzimuth,frange)

%%BeamPattern 
for idx_freq=1:length(frange)
    % Evaluate the beampattern and limit to max -40dB
    BPatternRFSB(idx_freq,:) = max(-40,mag2db(abs(squeeze(steerV(idx_freq,:,:))*Wrfsb(:,idx_freq))));
    BPatternMVDR(idx_freq,:) = max(-40,mag2db(abs(squeeze(steerV(idx_freq,:,:))*Wmvdr(:,idx_freq))));
end
figure(1)
subplot(1,2,1);
imagesc(angRangeAzimuth,frange(8:end),BPatternMVDR(8:end,:))
xlabel('DOA in degrees')
ylabel('Frequency in Hz')
title('Beampattern MVDR')
colorbar

subplot(1,2,2);
imagesc(angRangeAzimuth,frange(8:end),BPatternRFSB(8:end,:))
xlabel('DOA in degrees')
ylabel('Frequency in Hz')
title('Beampattern RFSB')
colorbar
%%