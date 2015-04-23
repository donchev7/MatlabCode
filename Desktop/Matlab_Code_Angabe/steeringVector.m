function steerVec = steeringVector(cfg,mic_position,frequencies)
cfg.fstep = 10; % frequency steps between 'sampling points' of set of
% lower and upper limit in order to create extended frequency range
% cfg.lfreq_lim = 100; %lower limit
% cfg.hfreq_lim = 200; %to create higher limit
% cfg.frange = 300:cfg.fstep:(cfg.fs/2 - cfg.hfreq_lim);
% cfg.frange = (cfg.frange(1)-cfg.lfreq_lim):cfg.fstep:(cfg.frange(end)+cfg.hfreq_lim);
cfg.frange = frequencies(:,1);
cfg.k_range = 2*pi*cfg.frange/cfg.c;
% 
% mic_position.x = [-7.5e-2 -3.375e-2 0 3.375e-2 7.5e-2];
% mic_position.y = [-6e-2 -0.75e-2 0 -0.75e-2 -6e-2];
% mic_position.z = [-4e-2 0 0 0 -4e-2]; 
cfg.mic_pos.x = mic_position.x;
cfg.mic_pos.y = mic_position.y;
cfg.mic_pos.z = mic_position.z;
steerVec = zeros(length(cfg.frange),cfg.nmic);
for idx_mic=1:cfg.nmic
    pos_mic = [cfg.mic_pos.x(idx_mic); cfg.mic_pos.y(idx_mic); cfg.mic_pos.z(idx_mic)];
    for idx_freq=1:length(cfg.frange)
        
%cfg.angular_resolution.azimuth = cfg.theta_vec;
%cfg.angular_resolution.elevation = repmat(cfg.look_elevation,size(cfg.angular_resolution.azimuth));

        kvec = - cfg.k_range(idx_freq) * [sind(cfg.look_elevation).*cosd(cfg.look_azimuth);...
        sind(cfg.look_elevation).*sind(cfg.look_azimuth); cosd(cfg.look_elevation)];
%create steering vector according to van Trees (2.28) for all
%microphone positions
        steerVec(idx_freq,idx_mic) = exp(-1i*kvec.'*pos_mic);
        %lam=cfg.c/cfg.frange(idx_freq);
        %steerVec(idx_freq,idx_mic)=steervec(pos_mic/lam,[cfg.look_azimuth;cfg.look_elevation]);
    end
end

end