function newW = interpolateFrequencies(frange_ext,cfg,flt)
    for idx_freq=1:length(frange_ext)
        if(idx_freq) >= 2
            k_new=find(cfg.frange <= frange_ext(idx_freq));
            C = setxor(k,k_new);
            for idx_mic=1:cfg.nmic
                newW(idx_mic,C) = flt.w.RFSB(idx_mic,idx_freq);
            end
            k=k_new;
        else
        k=find(cfg.frange <= frange_ext(idx_freq));
        for idx_mic=1:cfg.nmic
            newW(idx_mic,k) = flt.w.RFSB(idx_mic,idx_freq);
        end
        end
    end
end