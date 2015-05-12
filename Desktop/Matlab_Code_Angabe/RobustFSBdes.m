% RobustFSB Designs the optimum filter weights of a robust least-squares
% frequency-invariant polynomial (RLSFIP) filter-and-sum beamformer.
%

function [fir_imp_resp, cfg, steerVector, realWNG_dB] = ...
    RobustFSBdes(cfg)

    local = [];
    
    %----------------------------------------------------------------------
    % Initialization
    %----------------------------------------------------------------------
    %cfg.spacing = spacing; %spacing of microphones (or radius)  
    
    %----------------------------------------------------------------------
    % other parameters
  
    %----------------------------------------------------------------------
    % set path of set of prototype hrirs if required
    if strcmp(cfg.design, 'hrtf')
        cfg.path_hrirs = 'data/HRIR_NAO_LRC_az_0_355_16kHz.mat';
    end
    %----------------------------------------------------------------------
    local.WNG_lim = 10*log10(cfg.nmic); %limit of WNG
    if (cfg.wng_limit_db > local.WNG_lim)
        disp(' ')
        disp('White noise gain cannot be larger than 10*log10(N_transducer) corresponding to a delay-and-sum beamformer.')
        eval(['disp(''Reverting to maximum possible WNG: ' num2str(cfg.WNG_lim) 'dB'')'])
        disp(' ')
        cfg.wng_limit_db = local.WNG_lim;                
    end
    % Converts WNG into an absolute number    
    local.WNG = 10^(cfg.wng_limit_db/10);
    cfg.nmbr_look_dir = length(cfg.des_look_dir); %number of prototype look
    % directions (for simple FSB nmbr_look_dir =1)
    %----------------------------------------------------------------------
    % Initialisation of geometry dependent parameters        
    cfg = BF_Array_Geometry(cfg);
    
    %%
    %--------------------------------------------------------------
    % Initialisation of desired magnitude response
    % angular resolution (discretized angles) in degrees Todo: in which direction
    % be careful with changing this: I think cfg.desResponse_def is created for steps
    % of 5° (H.Barfuss, 13.11.14)
    local.P = 0;
    local.D_i = (cfg.des_look_dir.azimuth - 90)/90;          

    %local.angular_resolution.azimuth = cfg.angRange;
    %local.angular_resolution.elevation = repmat(cfg.des_look_dir.elevation,size(local.angular_resolution.azimuth));
    % define shape of desired magnitude response
    local.resp_choice = 'narrow';     % 'narrow', 'wide'
    switch local.resp_choice
        case 'narrow'
            local.desResponse_def = [linspace(0,0.7079,4) 0.99 1 0.99 fliplr(linspace(0,0.7079,4))];
        case 'wide'
            local.desResponse_def = [linspace(eps,1-eps,6) 1 linspace(1-eps,eps,6)];
    end
    %needed to determine starting and ending index (lower and upper
    %limit) of desired impulse response
    offset = floor(length(local.desResponse_def)/2);
    %contains the desired response for each prototype look
    %direction (column-wise)
    local.desResponse = zeros(length(cfg.angRange.azimuth),length(cfg.des_look_dir.azimuth));
    for idx_look_dir = 1:cfg.nmbr_look_dir 
        % skip_u and skip_l indicate if there was a problem with the
        % upper and lower limit of the desired response. If there
        % is a problem, one/both is > 0 and includes the difference
        % between maximum/minimum possible index and actual index
        skip_u = 0;
        skip_l = 0;
        %lower limit (starting index) of desired response
        lower_limit = find(cfg.angRange.azimuth == cfg.des_look_dir.azimuth(idx_look_dir)) - offset;
        if lower_limit < 1 % detects lower limit problem
            skip_l = abs(lower_limit-1);
            lower_limit = 1;
        end
        %upper limit (ending index) of desired response
        upper_limit = (lower_limit+length(local.desResponse_def)-1-skip_l);          
        if upper_limit > length(cfg.angRange.azimuth) % detects upper limit problem
            skip_u = upper_limit - length(cfg.angRange.azimuth);
            upper_limit = length(cfg.angRange.azimuth);
        end
        % computing indices
        ind_temp2 = lower_limit:upper_limit;      
        %check if lower or upper limit exceeds minimum/maximum
        %value, and store des desired reponse in the corresponding
        %column of cfg.desResponse
        if skip_l ~= 0
            local.desResponse(ind_temp2,idx_look_dir) = local.desResponse_def(skip_l+1:end);                             
        end
        if skip_u ~= 0
            local.desResponse(ind_temp2,idx_look_dir) = local.desResponse_def(1:end-skip_u);                             
        end
        if skip_l == 0 && skip_u == 0
              local.desResponse(ind_temp2,idx_look_dir) = local.desResponse_def;
        end
    end
    %--------------------------------------------------------------
    % initialisation of array response matrix G (see Equation (4) 
    % in [Mabande et al, 2009] at desired response angles 
    % [frequency,theta_resolution,N] = [#frequencies, #angles, #microphones]
    switch cfg.design
        case 'freefield'
            for idx_micPos=1:cfg.nmic
                %                 for cnt_LD = 1 : cfg.nmbr_look_dir %so far as I can see,
                %                 this for loop is not needed here (at the moment), maybe
                %                 needed for a polynomial beamformer...
                %_ext means that G_ext is created using the extended frequency
                %range (see RobustFSB.m)
                %x_mic: position vector of idx_micPos-th microphone
                mic_position = [cfg.mic_pos.x(idx_micPos); cfg.mic_pos.y(idx_micPos); cfg.mic_pos.z(idx_micPos)];
                %compute steering vectors and store them in cfg.G_ext
                for idx_frequency = 1:length(cfg.k_range)
                    %k_vec: wavevectors for current frequency i_frequency
                    %and current microphone with respect to all source positions
                    k_vec = - cfg.k_range(idx_frequency) * [sind(cfg.angRange.elevation).*cosd(cfg.angRange.azimuth); sind(cfg.angRange.elevation).*sind(cfg.angRange.azimuth); cosd(cfg.angRange.elevation)];
                    %save steering vectors in cfg.G_ext
                    % same as D(f,æ) eq. 30 Mic. Arrays tut
                    local.G_ext(idx_frequency,:,idx_micPos) = exp(-1i*k_vec.'*mic_position);
                end
                %                 end
            end
        case 'hrtf'
            % load hrirs
            load(cfg.path_hrirs);
            local.hrirs = imp_resp;
            %transform them into dft domain
            local.hrtfs = fft(imp_resp,cfg.fs,1);
            %save hrtfs in cfg.G_ext
            local.G_ext = permute(cfg.hrtfs(cfg.frange,cfg.idx_hrtfs,...
                cfg.angRange.azimuth/5+1),[1,3,2]);
    end
    %%
    
    %----------------------------------------------------------------------
    % Initialisation of look direction vector            
    local = InitLookDirVec(cfg,local);
    N_f_bins = length(cfg.k_range);
    %cfg.calG (\mathcal{G} in [Mabande et al, 2010]) includes the product 
    %G(\omega_p)*D_i of Eq. (5) for each (design-) frequency \omega_p
    %stacked in vertical direction for each prototype look direction
    %dimension: [Num_DesiredLookDirections*theta_resolution,Num_mics*(P+1)]
    %(see definition of \mathcal{G} in [Mabande et al, 2010] between Eq. (5) and (6))
    local.calG = [];
    %G_D is a temporary variable and will include the product G(\omega_p)*D_i 
    %from Eq. (5) of dimension: [freq,Num_DesiredLookDirections,theta_resolution,Num_mics*(P+1)]
    %this product is computed for every (design-) frequency
    G_D = NaN(N_f_bins, 1, length(local.desResponse), cfg.nmic);
    %for-loop over each (design-) frequency
    for idx_frequency = 1:length(cfg.frange)     
		tmp3 = [];    
        % calculation of G_D at frequency \omega_p of dimension [theta_resolution,Num_mics*(P+1)]
        for idx_look_dir = 1:cfg.nmbr_look_dir 
            %cfg.G_ext = steering vector
	        G_D(idx_frequency,idx_look_dir,:,:) = squeeze(local.G_ext(idx_frequency,:,:,idx_look_dir))*...
                squeeze(local.D(idx_look_dir,:,:));
        	% concatenating of matrices along first dimension ->
    	    % frequency-dependent G_D(\omega_p) are concatenated along vertical
	        % direction        
        	tmp3 = cat(1,tmp3,squeeze(G_D(idx_frequency,idx_look_dir,:,:)));         	 
        end
        %final result is stored in cfg.calG
        local.calG(idx_frequency,:,:) = tmp3;
    end
    % store all desired response (different for each desired prototype look
    % direction) in cfg.b_des (see Eq. (6) and definition before (6) in 
    % [Mabande et al, 2010])
    %
%     for cnt1 = 1:cfg.nmbr_look_dir %I think this for loop is not necessary
        local.b_des = local.desResponse(:);
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
    flt.w_opt = NaN(cfg.nmic,N_f_bins);
    %for loop over each (design-) frequency
    for idx_frequency=1:length(cfg.k_range)
        
        % ai_Di will include the product a^T_i(\omega_p)*D_i of dimension [1, N(P+1)]
        % on each row and is required for the distortionless constraint
        % -> Eq. (7) in [Mabande et al, 2010]        
        %
        % Thus, ai_Di is of dimension [I, N(P+1)] (I: number of prototype
        % directions)
        %
        % Note: ai_Di is also already included in cfg.calD at position
        % squeeze(cfg.calG(idx_frequency,cfg.angular_resolution == cfg.des_look_dir(i),:))
        ai_Di = [];
        for idx_look_dir = 1:cfg.nmbr_look_dir
            ai_Di(idx_look_dir,:) = squeeze(G_D(idx_frequency, idx_look_dir, cfg.angRange.azimuth == cfg.des_look_dir.azimuth(idx_look_dir),:));
        end
        
        % array responses at freq idx_frequency
        A = squeeze(local.calG(idx_frequency,:,:));
        
        %------------------------------------------------------------------
        % begin cvx optimization
        %------------------------------------------------------------------
        cvx_begin
        % pause for 1 sec
        cvx_quiet(1)     
        % complex variable declaration 
        variable w(cfg.nmic*(local.P+1),1) complex;          

        % minization of least-squares (LS) problem 
        % -> Eq. (6) in [Mabande et al, 2010]
        minimize ( norm(A*w - local.b_des, 2) );      
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
                norm(squeeze(local.D(idx_look_dir,:,:))*w,2) <= 1/sqrt(local.WNG);
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
    
%---------------------------------------------------------------------
%                       WNG calculations                              %
%---------------------------------------------------------------------
local.plot_LookDir = 1;
for idx_freq=1:length(cfg.frange) 
    % Steering Vector of dimension [#Microphones, 1]
    d = squeeze(G_D(idx_freq,local.plot_LookDir, cfg.angRange.azimuth == cfg.des_look_dir.azimuth(local.plot_LookDir),:)); 
    % theoretical WNG (computed with optimum filter coefficients
    % flt.w_opt obtained from optimization problem
    local.WNG(idx_freq) = 10*log10( (abs(flt.w_opt(:,idx_freq).'*d))^2/...
        ((squeeze(local.D(local.plot_LookDir,:,:))*flt.w_opt(:,idx_freq))'*(squeeze(local.D(local.plot_LookDir,:,:))*flt.w_opt(:,idx_freq))) );      
end 



%----------------------------------------------------------------------
%                           Set return values                         %
%----------------------------------------------------------------------
fir_imp_resp = flt.w_opt;
% Resulting WNG of designed beamformer in dB
realWNG_dB = local.WNG(:);
steerVector = local.G_ext;

end