%--------------------------------------------------------------------------
%                        Array Geometry                                   
%--------------------------------------------------------------------------
function cfg = BF_Array_Geometry(cfg) 
    
    switch (cfg.geometry)
        %------------------------------------------------------------------
        % define array geometry (geometry-dependent variables) for linear array
        % positions of microphones defined in cartesian coordinate systems
        % according to [Van Trees, Optimum Array Processing, Fig.2.1]
        % x: right(>0), left(<0)
        % y: forward(>0), backward(<0)
        % z: above(>0), below(<0)
        % origin of array (x,y,z)=(0,0,0) assumed to be center microphone
        %------------------------------------------------------------------
        case 1 % linear array
            if cfg.spacing == 0 %input non-uniform spacing manually
                switch cfg.design
                    case 'freefield'
                        switch cfg.N
                            case 3 %left, right, and center microphone (mics: 1,5, and 9)
                                %microphone positions always from left (mic 1) to
                                %right (mic 9)
                                cfg.mic_pos.z = mic_positions.z(1:3);                               
                                cfg.mic_pos.x = [-7.5e-2 0 7.5e-2];
                                cfg.mic_pos.y = [-6e-2 0 -6e-2];
                                cfg.mic_pos.z = [-4e-2 0 -4e-2];                                
                            case 5 %(mics: 1, 3, 5, 7, and 9)
                                cfg.mic_pos.x = [-7.5e-2 -3.375e-2 0 3.375e-2 7.5e-2];
                                cfg.mic_pos.y = [-6e-2 -0.75e-2 0 -0.75e-2 -6e-2];
                                cfg.mic_pos.z = [-4e-2 0 0 0 -4e-2];                                
                            case 7 %(mics: 1, 2, 4, 5, 6, 8, and 9)
                                cfg.mic_pos.x = [-7.5e-2 -4.5e-2 -2.25e-2 0 2.25e-2 4.5e-2 7.5e-2];
                                cfg.mic_pos.y = [-6e-2 -1e-2 -0.5e-2 0 -0.5e-2 -1e-2 -6e-2];
                                cfg.mic_pos.z = [-4e-2 0 0 0 0 0 -4e-2];                                
                            case 9 %(mics: 1, 2, 3, 4, 5, 6, 7, 8, and 9)
                                cfg.mic_pos.x = [-7.5e-2 -4.5e-2 -3.375e-2 -2.25e-2 0 2.25e-2 3.375e-2 4.5e-2 7.5e-2];
                                cfg.mic_pos.y = [-6e-2 -1e-2 -0.75e-2 -0.5e-2 0 -0.5e-2 -0.75e-2 -1e-2 -6e-2];
                                cfg.mic_pos.z = [-4e-2 0 0 0 0 0 0 0 -4e-2];                                
                        end
                    case 'hrtf'
                        switch cfg.N
                            case 3
                                cfg.idx_hrtfs = [1 5 9];
                            case 5
                                cfg.idx_hrtfs = [1 3 5 7 9];
                            case 7
                                cfg.idx_hrtfs = [1 2 4 5 6 8 9];
                            case 9
                                cfg.idx_hrtfs = (1:9);
                        end                         
                end
            else
                % initialisation of sensor positions
                cfg.mic_pos.x = linspace(-(cfg.N - 1)/2,(cfg.N - 1)/2,cfg.N)*cfg.spacing; 
                cfg.mic_pos.y = zeros(size(cfg.mic_pos.x));
                cfg.mic_pos.z = zeros(size(cfg.mic_pos.x));
            end

            %--------------------------------------------------------------
            % filterlength
            cfg.filterlength = 1023;                         
            %--------------------------------------------------------------
            % Desired prototype look directions (MAY HAVE TO BE ADAPTED FOR
            % POLYNOMIAL BEAMFORMING)
            %evtl. kann cfg.DesAng auch ganz entfernt werden!
            cfg.DesAng = cfg.des_look_dir;
            % range of delays required for polynomial beamforming (-> one
            % value for normal FSB)
            % NOTE: at the moment, elevation is always assumed to be fix.
            % Steering the beamformer is only done in azimuth direction!
            %das hier muss noch für den polynomial beamformer angepasst
            %werden denke ich!
            cfg.D_i = (cfg.DesAng.azimuth - 90)/90;          
            %--------------------------------------------------------------
            % Initialisation of desired magnitude response
            % angular resolution (discretized angles) in degrees Todo: in which direction
            % be careful with changing this: I think cfg.desResponse_def is created for steps
            % of 5° (H.Barfuss, 13.11.14)
            cfg.angular_resolution.azimuth = (0:5:180);
            cfg.angular_resolution.elevation = repmat(cfg.DesAng.elevation,size(cfg.angular_resolution.azimuth));
            % define shape of desired magnitude response
            cfg.resp_choice = 'narrow';     % 'narrow', 'wide'
            switch cfg.resp_choice
                case 'narrow'
                    cfg.desResponse_def = [linspace(0,0.7079,4) 0.99 1 0.99 fliplr(linspace(0,0.7079,4))];
                case 'wide'
                    cfg.desResponse_def = [linspace(eps,1-eps,6) 1 linspace(1-eps,eps,6)];
            end
            %needed to determine starting and ending index (lower and upper
            %limit) of desired impulse response
            offset = floor(length(cfg.desResponse_def)/2);
            %contains the desired response for each prototype look
            %direction (column-wise)
            cfg.desResponse = zeros(length(cfg.angular_resolution.azimuth),length(cfg.DesAng.azimuth));
            for idx_look_dir = 1:cfg.nmbr_look_dir 
                % skip_u and skip_l indicate if there was a problem with the
                % upper and lower limit of the desired response. If there
                % is a problem, one/both is > 0 and includes the difference
                % between maximum/minimum possible index and actual index
                skip_u = 0;
                skip_l = 0;
                %lower limit (starting index) of desired response
                lower_limit = find(cfg.angular_resolution.azimuth == cfg.DesAng.azimuth(idx_look_dir)) - offset;
                if lower_limit < 1 % detects lower limit problem
                    skip_l = abs(lower_limit-1);
                    lower_limit = 1;
                end
                %upper limit (ending index) of desired response
                upper_limit = (lower_limit+length(cfg.desResponse_def)-1-skip_l);          
                if upper_limit > length(cfg.angular_resolution.azimuth) % detects upper limit problem
                    skip_u = upper_limit - length(cfg.angular_resolution.azimuth);
                    upper_limit = length(cfg.angular_resolution.azimuth);
                end
                % computing indices
                ind_temp2 = lower_limit:upper_limit;      
                %check if lower or upper limit exceeds minimum/maximum
                %value, and store des desired reponse in the corresponding
                %column of cfg.desResponse
                if skip_l ~= 0
                    cfg.desResponse(ind_temp2,idx_look_dir) = cfg.desResponse_def(skip_l+1:end);                             
                end
                if skip_u ~= 0
                    cfg.desResponse(ind_temp2,idx_look_dir) = cfg.desResponse_def(1:end-skip_u);                             
                end
                if skip_l == 0 && skip_u == 0
                      cfg.desResponse(ind_temp2,idx_look_dir) = cfg.desResponse_def;
                end
            end
            %--------------------------------------------------------------
            % initialisation of array response matrix G (see Equation (4) 
            % in [Mabande et al, 2009] at desired response angles 
            % [frequency,theta_resolution,N] = [#frequencies, #angles, #microphones]
            switch cfg.design
                case 'freefield'
                    for idx_micPos=1:cfg.N
                        %                 for cnt_LD = 1 : cfg.nmbr_look_dir %so far as I can see,
                        %                 this for loop is not needed here (at the moment), maybe
                        %                 needed for a polynomial beamformer...
                        %_ext means that G_ext is created using the extended frequency
                        %range (see RobustFSB.m)
                        %x_mic: position vector of idx_micPos-th microphone
                        x_mic = [cfg.mic_pos.x(idx_micPos); cfg.mic_pos.y(idx_micPos); cfg.mic_pos.z(idx_micPos)];
                        %compute steering vectors and store them in cfg.G_ext
                        for idx_frequency = 1:length(cfg.k_range_ext)
                            %k_vec: wavevectors for current frequency i_frequency
                            %and current microphone with respect to all source positions
                            k_vec = - cfg.k_range_ext(idx_frequency) * [sind(cfg.angular_resolution.elevation).*cosd(cfg.angular_resolution.azimuth); sind(cfg.angular_resolution.elevation).*sind(cfg.angular_resolution.azimuth); cosd(cfg.angular_resolution.elevation)];
                            %save steering vectors in cfg.G_ext
                            cfg.G_ext(idx_frequency,:,idx_micPos) = exp(-1i*k_vec.'*x_mic);
                        end
                        %                 end
                    end
                case 'hrtf'
                    % load hrirs
                    load(cfg.path_hrirs);
                    cfg.hrirs = imp_resp;
                    %transform them into dft domain
                    cfg.hrtfs = fft(imp_resp,cfg.srate,1);
                    %save hrtfs in cfg.G_ext
                    cfg.G_ext = permute(cfg.hrtfs(cfg.frange_ext,cfg.idx_hrtfs,...
                        cfg.angular_resolution.azimuth/5+1),[1,3,2]);
            end
            %--------------------------------------------------------------
            
        %------------------------------------------------------------------
        % define array geometry (geometry-dependent variables) for circular array
        %------------------------------------------------------------------            
        case 2 
            error('Other Geometry not supported yet\n')
    end
end