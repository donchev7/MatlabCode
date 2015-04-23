% RobustFSB Designs the optimum filter weights of a robust least-squares
% frequency-invariant polynomial (RLSFIP) filter-and-sum beamformer.
%
% Inputs:
%                  N:   Number of sensors or actors
%            spacing:   spacing between micophones or radius of array [meters]
%                           0 => non-uniform spacing
%             WNG_dB:   Desired White Noise Gain in dB 
%                  P:   polynomial order = 0 for FSB
%         int_choice:   Interpolation choice 1=>polynomial interpolation; 2=> ??? interpolation  
%        norm_choice:   1 => L_2 norm; 2 => L_inf norm
%           geometry:   1=> linear; 2=> circular
%       des_look_dir:   desired look direction
%       design: 'freefield ' or 'hrtf'
%
% Outputs: TODO
%       fir_imp_resp:       Filter impulse response of approximated FIR
%                           filters
%       G_plot:             Matrix of steering vectors for look direction
%                           for which the beampattern has been plotted
%       resp:               Contains the magnitude of the beamformer response
%                           range is between 0 and 1
%       FrequencyAxis:      frequency bins used for the plots
%       AngularAxis:        angles used for the plots
%       realWNG_dB:         Resulting WNG of designed beamformer in dB
%
% Description:		
% This programm computes the optimum filter tabs of a broadband polynomial 
% beamformer. The goal is to minimize the sum of all squares subject to a 
% distortionless and WNG constraint. 
% CVX is used to solve the convex Second order Cone Program (SOCP).
%
% Edwin Mabande, Erlangen 01.11 

function [fir_imp_resp, G_plot, resp, faxis, angaxis, realWNG_dB] = ...
    RobustFSB(N,spacing,WNG_dB,P,int_choice,norm_choice,geometry, des_look_dir, design)

    if nargin ~= 9,
      error('9 arguments required');
    end  
    
    %----------------------------------------------------------------------
    % Initialization
    %----------------------------------------------------------------------
    cfg = [];
    % save arguments in the cfg struct
    cfg.N = N; %number of microphones
    cfg.spacing = spacing; %spacing of microphones (or radius)
    cfg.P = P; %order of polynomial beamformer, if 0 -> simple FSB
    cfg.int_choice = int_choice; %choice of interpolation, not used in this code!
    cfg.norm_choice = norm_choice; %choice of norm to be minimized    
    cfg.geometry = geometry; %parameter indicating array geometry (linear or spherical)
    cfg.des_look_dir = des_look_dir; %desired look direction
    cfg.design = design; %choice of freefield-based or hrtf-based design
    %----------------------------------------------------------------------
    % other parameters
    cfg.c = 342; %speed of sound
    cfg.srate = 16000; %sampling rate
    %----------------------------------------------------------------------
    % set path of set of prototype hrirs if required
    if strcmp(cfg.design, 'hrtf')
        cfg.path_hrirs = 'data/HRIR_NAO_LRC_az_0_355_16kHz.mat';
    end
    %----------------------------------------------------------------------
    % specify beamforming frequency range
    cfg.fstep = 100; % frequency steps between 'sampling points' of set of 
                     % discrete frequencies used to design the filter weights
    %
    % lower and upper limit in order to create extended frequency range
    cfg.lfreq_lim = 100; %lower limit
    cfg.hfreq_lim = 200; %to create higher limit
    %
    % set of frequencies used to design the filter weights
    %Note: Upper Bound + cfg.hfreq_lim <= cfg.srate, otherwise fir2 (below)
    %won't work
    cfg.frange = 300:cfg.fstep:(cfg.srate/2 - cfg.hfreq_lim);
    %
    % extended frequency range for computations: extends beamforming frequency range 
    % computation in order to ensure a flat response across entire frequency range!
    cfg.frange_ext = (cfg.frange(1)-cfg.lfreq_lim):cfg.fstep:(cfg.frange(end)+cfg.hfreq_lim); 
    %
    cfg.k_range_ext = 2*pi*cfg.frange_ext/cfg.c; % wavenumber (corresponding to extended computation frequency range)
    cfg.k_range = 2*pi*cfg.frange/cfg.c; % wavenumber (corresponding to beamformer frequency range)
    %        
    cfg.WNG_dB = WNG_dB;
    cfg.WNG_lim = 10*log10(N); %limit of WNG
    if (cfg.WNG_dB > cfg.WNG_lim)
        disp(' ')
        disp('White noise gain cannot be larger than 10*log10(N_transducer) corresponding to a delay-and-sum beamformer.')
        eval(['disp(''Reverting to maximum possible WNG: ' num2str(cfg.WNG_lim) 'dB'')'])
        disp(' ')
        cfg.WNG_dB = cfg.WNG_lim;                
    end
    % Converts WNG into an absolute number    
    cfg.WNG = 10^(cfg.WNG_dB/10);
    %    
    % Possibly NEEDS TO BE ADAPTED FOR POLYNOMIAL BEAMFORMER
    cfg.nmbr_look_dir = length(cfg.des_look_dir); %number of prototype look
    % directions (for simple FSB -> 1)
    %----------------------------------------------------------------------
    % Initialisation of geometry dependent parameters        
    cfg = BF_Array_Geometry(cfg);
    %----------------------------------------------------------------------
    % Initialisation of look direction vector            
    cfg = InitLookDirVec(cfg);
    %----------------------------------------------------------------------
    % Initialising of matrices for the optimization problem
    %
    % ToDo: I think this needs to be adapted (regarding the concatenation)
    % for several desired look directions -> for polynomial beamforming
    % H.Barfuss, 10.11.14
    %
    N_f_bins = length(cfg.k_range_ext); %number of (design-) frequencies
    %cfg.calG (\mathcal{G} in [Mabande et al, 2010]) includes the product 
    %G(\omega_p)*D_i of Eq. (5) for each (design-) frequency \omega_p
    %stacked in vertical direction for each prototype look direction
    %dimension: [Num_DesiredLookDirections*theta_resolution,Num_mics*(P+1)]
    %(see definition of \mathcal{G} in [Mabande et al, 2010] between Eq. (5) and (6))
    cfg.calG = [];
    %G_D is a temporary variable and will include the product G(\omega_p)*D_i 
    %from Eq. (5) of dimension: [freq,Num_DesiredLookDirections,theta_resolution,Num_mics*(P+1)]
    %this product is computed for every (design-) frequency
    G_D = NaN(N_f_bins, 1, length(cfg.desResponse), N);
    %for-loop over each (design-) frequency
    for idx_frequency = 1:length(cfg.frange_ext)     
		tmp3 = [];    
        % calculation of G_D at frequency \omega_p of dimension [theta_resolution,Num_mics*(P+1)]
        for idx_look_dir = 1:cfg.nmbr_look_dir 
            %cfg.G_ext = steering vector
	        G_D(idx_frequency,idx_look_dir,:,:) = squeeze(cfg.G_ext(idx_frequency,:,:,idx_look_dir))*...
                squeeze(cfg.D(idx_look_dir,:,:));
        	% concatenating of matrices along first dimension ->
    	    % frequency-dependent G_D(\omega_p) are concatenated along vertical
	        % direction        
        	tmp3 = cat(1,tmp3,squeeze(G_D(idx_frequency,idx_look_dir,:,:)));         	 
		end
        %final result is stored in cfg.calG
        cfg.calG(idx_frequency,:,:) = tmp3;
    end
    % store all desired response (different for each desired prototype look
    % direction) in cfg.b_des (see Eq. (6) and definition before (6) in 
    % [Mabande et al, 2010])
    %
%     for cnt1 = 1:cfg.nmbr_look_dir %I think this for loop is not necessary
        cfg.b_des = cfg.desResponse(:);
%     end
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    %                     CVX based Optimization                          %
    %----------------------------------------------------------------------    
    % This block of code performs the constraint optimization proposed in 
    % [Mabande et al, 2010] to obtain the optimum filter weighs in the
    % design domain
    % The CVX toolbox is used to solve the (convex) optimization problem
    %----------------------------------------------------------------------
%     % Use a waitbar to illustrate design progress
%     h_wb = waitbar(0, 'Beamformer is being designed...');
    %----------------------------------------------------------------------
    % variable for resulting optimum filter weights
    % This will include the optimum filter coefficients for the extended
    % design frequency range cfg.frange_ext
    flt.w_opt = NaN(cfg.N,N_f_bins);
    %for loop over each (design-) frequency
    for idx_frequency=1:length(cfg.k_range_ext)
        
        % ai_Di will include the product a^T_i(\omega_p)*D_i of dimension [1, N(P+1)]
        % on each row and is required for the distortionless constraint
        % -> Eq. (7) in [Mabande et al, 2010]        
        %
        % Thus, ai_Di is of dimension [I, N(P+1)] (I: number of prototype
        % directions)
        %
        % Note: ai_Di is also already included in cfg.calD at position
        % squeeze(cfg.calG(idx_frequency,cfg.angular_resolution == cfg.DesAng(i),:))
        ai_Di = [];
        for idx_look_dir = 1:cfg.nmbr_look_dir
            ai_Di(idx_look_dir,:) = squeeze(G_D(idx_frequency, idx_look_dir, cfg.angular_resolution.azimuth == cfg.DesAng.azimuth(idx_look_dir),:));
        end
        
        % array responses at freq idx_frequency
        A = squeeze(cfg.calG(idx_frequency,:,:));
        
        %------------------------------------------------------------------
        % begin cvx optimization
        %------------------------------------------------------------------
        cvx_begin
        % pause for 1 sec
        cvx_quiet(1)     
        % complex variable declaration 
        variable w(cfg.N*(cfg.P+1),1) complex;          

        % minization of least-squares (LS) problem 
        % -> Eq. (6) in [Mabande et al, 2010]
        minimize ( norm(A*w - cfg.b_des, 2) );      
        % constraints of LS problem
        subject to
            % for each prototype look direction one distortionless constraints 
            % -> Eq. (7) in [Mabande et al, 2010]
            for idx_look_dir = 1:cfg.nmbr_look_dir         
                imag(w.'*ai_Di(idx_look_dir,:).') == 0;                       
                real(w.'*ai_Di(idx_look_dir,:).') == 1; 
            end
            
%             % For a response <= 1
%             max(abs(A*w)) <= 1.0;   
            
            % for each prototype look direction one WNG constraint
            % -> Eq. (8) in [Mabande et al, 2010]
            for idx_look_dir = 1:cfg.nmbr_look_dir         
                norm(squeeze(cfg.D(idx_look_dir,:,:))*w,2) <= 1/sqrt(cfg.WNG);
            end

        cvx_end
        %------------------------------------------------------------------
        %end of cvx optimization           
        %------------------------------------------------------------------        
%         cvx_optval 

        % Concatenate optimum filter coefficiens of current frequency
        % flt.w_opt contains the optimum filter weights of one filter at each row
        flt.w_opt(:,idx_frequency) = w;                     
        
%         % Update waitbar
%         waitbar(idx_frequency/N_f_bins, h_wb)
    end
%     % Close waitbar
%     close(h_wb);
    

    %----------------------------------------------------------------------
    %                             FIR Approximation                       %
    %----------------------------------------------------------------------    
    % filter coefficient computation -> FIR approxmation of optimum filters
    % saved in flt.w_opt
    %----------------------------------------------------------------------
    % w_opt_wholeFreqRange will contain the optimum filter weights to be
    % approximated by fir2
    %
    % At frequencies that are not included in the extended set of design 
    % frequencies cfg.frange_ext, the eps value is used
    w_opt_wholeFreqRange = ones((cfg.P+1)*cfg.N,(cfg.srate/2)/cfg.fstep+1)*eps;
    % F is a vector of frequency points in the range from 0 to 1, where 1 
    % corresponds to the Nyquist frequency (half the sampling frequency)
    % required for fir2 function
    F = (0:cfg.fstep:cfg.srate/2)/(cfg.srate/2);
    % write optimum filter weights of each filter to variable w_opt_wholeFreqRange
    % cfg.frange_ext/cfg.fstep+1 makes sure that at entries of w_opt_wholeFreqRange 
    % corresponding to frequencies that are not included in the extended set of
    % design frequencies cfg.frange_ext (e.g. 0Hz) no optimum filger weight is stored
    % This is necessesary because fir2 requires information about the whole
    % frequency range between 0Hz ... fs/2Hz
    for cnt = 1:cfg.N*(cfg.P+1)
        w_opt_wholeFreqRange(cnt,cfg.frange_ext/cfg.fstep+1) = ...
            conj(flt.w_opt(cnt,:)); %Why conj() ? weiß wohl nur Edwin :-)
    end
    % tmp_firfilt: contains approximated filter coefficients for each
    % filter on one row
    % -> dimension: [N*(P+1),filt_len]
    tmp_firfilt = zeros((cfg.P+1)*cfg.N,cfg.filterlength);
    for idx_filters=1:cfg.N*(cfg.P+1)
        tmp_firfilt(idx_filters,:) = fir2(cfg.filterlength-1,F,...
            squeeze(w_opt_wholeFreqRange(idx_filters,:)));              
    end

    % ordering coefficients in P+1 filter-and-sum subblocks of dimension 
    % [Filter length, #Microphones]. cfg.PBF_firfilt is used later on in
    % BF_Plot_BP
    % PBF_firfilt: (P+1) x filterlength x Number of mics
    cfg.PBF_firfilt = zeros((cfg.P+1),cfg.filterlength,cfg.N);% [P+1,filt_len, N]
    for idx_mic = 1:cfg.N
        for idx_P = 1:cfg.P+1
            cfg.PBF_firfilt(idx_P,:,idx_mic) = tmp_firfilt((idx_mic-1)*(cfg.P+1)+idx_P,:);      
        end
    end
    %----------------------------------------------------------------------
    %Choose look direction chosen from prototype look directions for which
    %the beampattern shall be plotted    
    cfg.plot_LookDir = 1;                    
    %----------------------------------------------------------------------


    %----------------------------------------------------------------------
    %                                 WNG                                 %
    %----------------------------------------------------------------------
    % In this section, the White Noise Gain (WNG) is computed for the
    % original design frequency range
    % First it is computed for cfg.frange_ext
    % Second, the resulting WNG is limited to cfg.frange 
    %----------------------------------------------------------------------
	% computing the frequenc responses with obtained fiter coefficients 
    % (of approximated FIR filters) 
    filt_real= zeros(cfg.N*(cfg.P+1),size(flt.w_opt,2));
    for idx_filter=1:cfg.N*(cfg.P+1)
        %compute frequency response of approximated FIR filters
        H  = freqz( tmp_firfilt(idx_filter,:), 1, cfg.srate/2+1 );
        %store compute frequency response of approximated FIR filters in
        %filt_real, one filter in each row
        filt_real(idx_filter,:) = H(cfg.frange_ext);
    end
    %----------------------------------------------------------------------
    %Compute theoretical (-> WNG_theoretical) and real (-> WNG_real) WNGs
    % in the 'original' frequency range, not the extended frequency range
    %For the computation of the WNG, see e.g., Eq. (8) in [Mabande et al,
    %2010]
    WNG_theoretical_tmp = NaN(size(cfg.frange_ext));
    WNG_real_tmp = NaN(size(cfg.frange_ext));
    % for each 
    for idx_freq_ext=1:length(cfg.frange_ext) 
        % Steering Vector of dimension [#Microphones, 1]
        d = squeeze(G_D(idx_freq_ext,cfg.plot_LookDir, cfg.angular_resolution.azimuth == cfg.DesAng.azimuth(cfg.plot_LookDir),:)); 
        % theoretical WNG (computed with optimum filter coefficients
        % flt.w_opt obtained from optimization problem
        WNG_theoretical_tmp(idx_freq_ext) = 10*log10( (abs(flt.w_opt(:,idx_freq_ext).'*d))^2/...
            ((squeeze(cfg.D(cfg.plot_LookDir,:,:))*flt.w_opt(:,idx_freq_ext))'*(squeeze(cfg.D(cfg.plot_LookDir,:,:))*flt.w_opt(:,idx_freq_ext))) );
        % real/actual WNG (computed from approximated fir filter
        % coefficients in filt_real)
        WNG_real_tmp(idx_freq_ext) = 10*log10( (abs(filt_real(:,idx_freq_ext).'*d))^2/...
            ((squeeze(cfg.D(cfg.plot_LookDir,:,:))*filt_real(:,idx_freq_ext))'*(squeeze(cfg.D(cfg.plot_LookDir,:,:))*filt_real(:,idx_freq_ext)))  );
    end 
    %save only WNG in the 'original' frequency range (cfg.frange), not in
    %the extended frequency range
    cfg.WNG_theoretical = WNG_theoretical_tmp(cfg.lfreq_lim/cfg.fstep+1:end - cfg.hfreq_lim/cfg.fstep);
    cfg.WNG_real = WNG_real_tmp(cfg.lfreq_lim/cfg.fstep+1:end - cfg.hfreq_lim/cfg.fstep);
    %----------------------------------------------------------------------

    
    %----------------------------------------------------------------------
    %                             Beampatterns                            %
    %----------------------------------------------------------------------
    % In this section the resulting beampattern and WNG is plotted
    %----------------------------------------------------------------------
    % Illustrate results
    cfg = BF_Plot_BP(cfg);
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    %                           Set return values                         %
    %----------------------------------------------------------------------
    % Filter impulse response of approximated FIR filters
    fir_imp_resp = tmp_firfilt.';
    %Matrix of steering vectors for look direction for which the beampattern has been plotted
    G_plot = cfg.G_plot;
    %Magnitude of the beamformer response
    resp = cfg.BF_response_abs;        
    % Frequency bins used for the plots
    faxis = cfg.frange(:);    
    % Angles used for the plots
    angaxis = cfg.DOA_degs_plot(:);
    % Resulting WNG of designed beamformer in dB
    realWNG_dB = cfg.WNG_real(:);
end



