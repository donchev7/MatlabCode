function Ymwf = postBeamformerProcessing(X,Thi_Y,cfg,GammaIso,d,Wmvdr)

%Wmvdr = myMVDR2(squeeze(GammaIso(idx_freq,:,:)),d(:,idx_freq),cfg.beamFormer,cfg.nmic);
%Wmwf as defined in Kuklasinski et. al EQ 5b 
Wmwf = postFilter(Wmvdr,GammaIso,d,Thi_Y,cfg.nmic);
% Ymvdr(idx_freq,:) = squeeze(X(idx_freq,:,:)).'*Wmvdr(:,idx_freq);
% Ymwf(idx_freq,:) = squeeze(X(idx_freq,:,:)).'*Wmvdr(:,idx_freq)*Wmwf;
% Ydsb(idx_freq,:) = squeeze(X(idx_freq,:,:)).'*(d(:,idx_freq)/cfg.nmic);

end