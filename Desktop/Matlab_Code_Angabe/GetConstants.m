function d_hrtf_tmp = GetConstants(cfg)

%Calculate the sound field as in [1] Kuklasinski



%%Evaluating gammaIso the cylindrical sound field from HRIR measurments
% Also evaluating the steering vector d from the RIR measurments as
% mentioned in Kuklasinski paper
HRIR=load('./Data/RIRs/HRIR_NAO_LRC_az_0_355_16kHz_norm.mat');
%L=length(HRIR.imp_resp(:,1,1));
[Nx_orig,K_signal,Ang] = size(HRIR.imp_resp);
n_zeros = ceil(Nx_orig/cfg.N)*cfg.N - Nx_orig;
for idx_mic=1:cfg.nmic
    %RIR_Padded = [flt.h(1:1058,cfg.look_azimuth/5+1,idx_mic); zeros(padding_zeros,1)];
    %d(idx_mic,:) = fft(RIR_Padded,cfg.N);
    for idx_ang=1:size(HRIR.imp_resp,3)
        HRIR_Padded = [HRIR.imp_resp(:,cfg.mic_ch(idx_mic),idx_ang); zeros(n_zeros,1)];
        d_hrtf_tmp(:,idx_ang,idx_mic) = fft(HRIR_Padded,cfg.N,1);
    end
end

d_hrtf_tmp = d_hrtf_tmp(1:cfg.N/2+1,:,:)/Nx_orig;
% for idx_ang=1:size(HRIR.imp_resp,3)
%     d_hrtf(:,:,idx_ang) = cov(d_hrtf_tmp(:,:,idx_ang));
% end
% gammaIsoSum = zeros(cfg.nmic,cfg.nmic);
% for idx_ang=1:size(HRIR.imp_resp,3)
%      gammaIsoSum = gammaIsoSum+d_hrtf(:,:,idx_ang)*d_hrtf(:,:,idx_ang)';
% end
% gammaIso=gammaIsoSum./size(HRIR.imp_resp,3);

%% End of evalution of gammaIso and HRTF steering vector d_hrtf_tmp
end