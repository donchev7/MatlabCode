function steerVec = steeringVector(cfg)

steerVec = zeros(length(cfg.frange),length(cfg.angRange.elevation),cfg.nmic);
for idx_mic=1:cfg.nmic
    pos_mic = [cfg.mic_pos.x(idx_mic); cfg.mic_pos.y(idx_mic); cfg.mic_pos.z(idx_mic)];
    for idx_freq=1:length(cfg.frange)
        kvec = - cfg.k_range(idx_freq) * [sind(cfg.angRange.elevation).*cosd(cfg.angRange.azimuth);...
        sind(cfg.angRange.elevation).*sind(cfg.angRange.azimuth); cosd(cfg.angRange.elevation)];

        steerVec(idx_freq,:,idx_mic) = exp(-1i*kvec.'*pos_mic);
    end
end

end