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
                        switch cfg.nmic
                            case 3 %left, right, and center microphone (mics: 1,5, and 9)
                                %microphone positions always from left (mic 1) to
                                %right (mic 9)                             
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
                        switch cfg.nmic
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
            else %uniform spacing automatically
                % initialisation of sensor positions
                cfg.mic_pos.x = linspace(-(cfg.nmic - 1)/2,(cfg.nmic - 1)/2,cfg.nmic)*cfg.spacing; 
                cfg.mic_pos.y = zeros(size(cfg.mic_pos.x));
                cfg.mic_pos.z = zeros(size(cfg.mic_pos.x));
            end
                       

        case 2
            disp('not implemented yet');
    end
end