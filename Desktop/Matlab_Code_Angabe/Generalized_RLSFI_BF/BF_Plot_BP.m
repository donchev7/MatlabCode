%--------------------------------------------------------------------------
%                    Plotting PBF Array Beampattern                       %
%--------------------------------------------------------------------------
% This function plots the resulting beampattern of the polynomial
% beamformer as well as the corresponding WNG
%--------------------------------------------------------------------------
function cfg = BF_Plot_BP(cfg)
    %----------------------------------------------------------------------
    % filt_real = frequency response of approximated FIR filters
    % dimension: (P+1) x #Microphones x frequency bins
    filt_real= [];                                           % 
    for i_PBF_firfilt = 1:cfg.P+1
        for i_mics=1:cfg.N
            % compute frequency response of approximated FIR filter
            % cfg.PBF_firfilt contains filter coefficients of P+1 filter-and-sum subblocks
            % of dimension [Filter length, #Microphones]            
            % -> dimension of cfg.PBF_firfilt: [(P+1) x filter length x
            % #Microphones]
            H  = freqz(squeeze(cfg.PBF_firfilt(i_PBF_firfilt,:,i_mics)), 1, cfg.srate/2+1);    
            %save frequency respones in filt_real
            %here H(cfg.frange), since we only want to plot the beampattern
            % in the 'original frequency range', not the extended range
            filt_real(i_PBF_firfilt,i_mics,:) = H(cfg.frange); 
        end
    end
    %----------------------------------------------------------------------
   
    % ---------------------------------------------------------------------
    % create steering vectors for plotting the beampatterns
    for i_mics=1:cfg.N
        switch cfg.geometry
            case 1 % linear array
                %cfg.DOA_degs_plot: DoAs used to plot the beampattern in
                %degree
                cfg.DOA_degs_plot.azimuth = (0:5:180);
                cfg.DOA_degs_plot.elevation = repmat(cfg.DesAng.elevation,size(cfg.DOA_degs_plot.azimuth));
				switch cfg.design
                    case 'freefield'
                        %x_mic: position vector of idx_micPos-th microphone
                        x_mic = [cfg.mic_pos.x(i_mics); cfg.mic_pos.y(i_mics); cfg.mic_pos.z(i_mics)];
                        % G_plot: stores all steering vectors for current array geometry
                        % created using the angular resolution cfg.DOA_degs_plot used
                        % for plotting the beampattern
                        % dimension: [frequency bins x #Angles for plotting x #Microphones]
                        for idx_frequency = 1:length(cfg.k_range)
                            %k_vec: wavevectors for current frequency i_frequency
                            %and current microphone with respect to all source positions
                            k_vec = - cfg.k_range(idx_frequency) * [sind(cfg.DOA_degs_plot.elevation).*cosd(cfg.DOA_degs_plot.azimuth); sind(cfg.DOA_degs_plot.elevation).*sind(cfg.DOA_degs_plot.azimuth); cosd(cfg.DOA_degs_plot.elevation)];
                            %%save steering vectors
                            G_plot(idx_frequency,:,i_mics) = exp(-1i*k_vec.'*x_mic);
                        end
                    case 'hrtf'
                        % create matrix of hrtfs used to plot the beampattern (for
                        % original frequency range cfg.frange!)
                        G_plot = permute(cfg.hrtfs(cfg.frange,cfg.idx_hrtfs,...
                            cfg.DOA_degs_plot.azimuth/5+1),[1,3,2]);
                end
            case 2 % circular array
                %cfg.DOA_degs_plot: DoAs used to plot the beampattern in radians
                cfg.DOA_degs_plot.azimuth = (0:1:360);
                cfg.DOA_degs_plot.elevation = repmat(cfg.DesAng.elevation,size(cfg.DOA_degs_plot.azimuth));
                %%TODO noch nicht geprüft
                for idx_frequency = 1:length(cfg.k_range)
                    %k_vec: wavevectors for current frequency i_frequency
                    %and current microphone with respect to all source positions
                    k_vec = - cfg.k_range(idx_frequency) * [sind(cfg.DOA_degs_plot.elevation).*cosd(cfg.DOA_degs_plot.azimuth); sind(cfg.DOA_degs_plot.elevation).*sind(cfg.DOA_degs_plot.azimuth); cosd(cfg.DOA_degs_plot.elevation)];
                    %%save steering vectors
                    G_plot(idx_frequency,:,i_mics) = exp(-1i*k_vec.'*x_mic);
                end  
        end
    end
    %
    cfg.G_plot = G_plot;
    % ---------------------------------------------------------------------
    

    % ---------------------------------------------------------------------
    % Create beamformer response according to (1) in [Mabande et al, 2010]
    %
    % k_range: wavenumbers in defined range (e.g., cfg.frange = 300 - 3400 Hz);
    % BF_response will contain the beamformer response of the polynopmial
    % beamformer with respect to the chosen look direction for plotting indexed
    % by fg.plot_LookDir (which is set in RobustFSB.m)
    % dimension -> [frequency bins, #Angles for plotting]
    BF_response=zeros(length(cfg.k_range),length(cfg.DOA_degs_plot.azimuth));
    % for each FSU and microphone (P+1: number of filter-and-sum units (FSUs),
    % N: Number of microphones)
    for i_FSUs = 1:cfg.P+1
        tmp = [];
        for i_mics=1:cfg.N
            % d_lookDir: Includes includes all P+1 interpolation factors
            % corresponding to the chose plot-look direction
            %
            % cfg.plot_LookDir: index of a prototype look direction out of
            % array cfg.DesAng
            % cfg.d: contains the D_i.^p and is of size I x P+1, created in
            % InitLookDirVec.m
            d_lookDir = cfg.d(cfg.plot_LookDir,:);
            % temp includes the frequency responses of the approximated FIR
            % filters of the current FSU weighted with the corresponding
            % interpolation factor
            % dimension: [freq,length(cfg.DOA_degs_plot),N]
            tmp(:,:,i_mics)=repmat(squeeze(filt_real(i_FSUs,i_mics,:))*d_lookDir(i_FSUs),1,length(cfg.DOA_degs_plot.azimuth));
        end
        %sum(G_plot.*tmp,3) creates the beamformer response of the current
        %FSU, weighted with the corresponding interpolation factor
        %the sum over all FSUs (outer loop) creates the overall beamformer
        %response according to (1) in [Mabande et al, 2010]
        BF_response = BF_response + sum(G_plot.*tmp,3);
    end
    % save created beamformer response
    cfg.BF_response = BF_response;
    % ---------------------------------------------------------------------


    % ---------------------------------------------------------------------
    %Plot beampattern, response in look direction and WNGs
    % ---------------------------------------------------------------------
    % save look direction for which beampattern should be plotted
    des_angle = cfg.des_look_dir.azimuth(cfg.plot_LookDir);
    % take absolute value of beamformer response and normalize to maximum
    % value
    BF_response_abs=abs(BF_response);
    BF_response_abs=BF_response_abs/max(max(BF_response_abs));
    %compute beampattern
    BPattern=20*log10(BF_response_abs);
    %limit lower values to -40dB
    a=find(BPattern<-40);
    BPattern(a)=-40;
    %save beampattern in cfg structure
    cfg.BPattern = BPattern;
    cfg.BF_response_abs = BF_response_abs;
    
    if 0
        % ---------------------------------------------------------------------
        %plot beampattern (created from approximated FIR filters)
        figure
        subplot(2,2,[1 3]);
        imagesc(cfg.DOA_degs_plot.azimuth,cfg.frange,BPattern)
        xlabel('DOA in degrees')
        ylabel('Frequency in Hz')
        title('Beampattern (FIR approx)')
        colorbar
        % ---------------------------------------------------------------------
        %compute beamformer response in look direction
        BF_lookDir_abs=BF_response_abs(:,find(cfg.DOA_degs_plot.azimuth == round(des_angle)));
        BF_lookDir_abs_log = 20*log10(BF_lookDir_abs);
        %limit from below by -20dB
        a=find(BF_lookDir_abs_log<-20);
        BF_lookDir_abs_log(a)=-20;
        %save desired response in look direction
        cfg.BF_lookDir_abs_log = BF_lookDir_abs_log;
        % ---------------------------------------------------------------------
        %plot response in look direction
        subplot(2,2,2);
        semilogx(cfg.frange,BF_lookDir_abs_log,'LineWidth',2);
        axis([10^(log10(cfg.frange(1))) 10^(log10(cfg.frange(end))) -20 0.1]);
        grid on; xlabel('Frequency in Hz'); ylabel('Response in dB')
        title(['Response in look-direction \vartheta = ', num2str(des_angle)])    
    %     set(gca, 'xtick', [400:100:1000, 2000, 3000]);
    %     set(gca, 'xticklabel', {'' '500' '' '' '' '' '1000' '' '3000'});
        % ---------------------------------------------------------------------
        %plot theoretical wng, real wng, and difference between theoretical and
        %real wng
        subplot(2,2,4)
        % plot(cfg.frange,cfg.WNG_theoretical,'b')
        semilogx(cfg.frange,cfg.WNG_real,'LineWidth',2); hold on;
        semilogx(cfg.frange,cfg.WNG_theoretical,'--r','LineWidth',2);
        semilogx(cfg.frange,cfg.WNG_theoretical-cfg.WNG_real,'-.k','LineWidth',2);
        legend('real WNG','theoretical WNG', 'difference')
        grid on
        xlim([10^(log10(cfg.frange(1))) 10^(log10(cfg.frange(end)))]);
        xlabel ('Frequency in Hz')
        ylabel ('White Noise Gain in dB')
        ylim([min(cfg.WNG_real)-3, max( max(cfg.WNG_real), min(cfg.WNG_real)+10)])
        % legend('Expected White Noise Gain','White Noise Gain realized by the FIR-filter',4)
        hold on
    %     set(gca, 'xtick', [400:100:1000, 2000, 3000]);
    %     set(gca, 'xticklabel', {'' '500' '' '' '' '' '1000' '' '3000'});
    end
end