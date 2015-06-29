function [d_hrtf, d] = GetConstants(cfg,flt)

%Calculate the sound field as in [1] Kuklasinski



%%Evaluating gammaIso the cylindrical sound field from HRIR measurments
% Also evaluating the steering vector d from the RIR measurments as
% mentioned in Kuklasinski paper
HRIR=load('./Data/RIRs/HRIR_NAO_LRC_az_0_355_16kHz_norm.mat');
%HRIR.imp_resp = HRIR.imp_resp/max(max(abs(HRIR.imp_resp(:,cfg.mic_ch,cfg.position(1))))); 
%RIR=load(cfg.path_imp_resp);
%L=length(HRIR.imp_resp(:,1,1));
%[Nx_orig,K_signal,Ang] = size(HRIR.imp_resp);
%n_zeros = ceil(Nx_orig/cfg.N)*cfg.N - Nx_orig;
[k, l] = max(flt.h(:,cfg.position(1),cfg.ref));
for idx_mic=1:cfg.nmic
    for idx_ang=1:size(HRIR.imp_resp,3)
        %HRIR_Padded = [HRIR.imp_resp(:,cfg.mic_ch(idx_mic),idx_ang); zeros(n_zeros,1)];
        d_hrtf_tmp(:,idx_ang,idx_mic) = fft(HRIR.imp_resp(:,cfg.mic_ch(idx_mic),idx_ang),cfg.N,1);
        if(strcmp(cfg.RIRcond,'NAO600_4m'))
            d(:,idx_ang,idx_mic) = fft([flt.h(l-20:l+80,idx_ang,idx_mic); zeros(80,1)],cfg.N,1);
        elseif(strcmp(cfg.RIRcond,'NAO600_2m'))
            d(:,idx_ang,idx_mic) = fft([flt.h(l-40:l+140,idx_ang,idx_mic); zeros(80,1)],cfg.N,1);
        else
            d(:,idx_ang,idx_mic) = fft([flt.h(l-40:l+200,idx_ang,idx_mic); zeros(80,1)],cfg.N,1);
        end
        %d(:,idx_ang,idx_mic) = fft(flt.h(l-20:l+320,idx_ang,idx_mic),cfg.N,1);
    end
end

d_hrtf = d_hrtf_tmp(1:cfg.N/2+1,:,:);
d = d(1:cfg.N/2+1,:,:);
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