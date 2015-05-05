function [FilterImpulseResponse] = main_robust_FSBORG(cfg, look_azimuth, look_elevation, design, wng_limit_db)
    % main_robust_FSB     This is the main file to design a robust least-squares
    % frequency-invariant (RLSFI) filter-and-sum beamformer.
    % In this file the parameters for the function RobustFSB(...) are set and 
    % RobustFSB(...) is called to design the beamformer weights.
    % Tis code is mainly based on the following two papers:
    %   * [Mabande et al,2009]: Design of Robust Superdirective Beamformers as
    %                           a Convex Optimization Problem, E. Mabande, 
    %                           A. Schad, and W. Kellermann
    %   * [Mabande et al,2010]: Design of Robust Polynomial Beamformers as a
    %                           Convex Optimization Problem, E. Mabande and W.
    %                           Kellermann

    %--------------------------------------------------------------------------
    % Initialization
    %--------------------------------------------------------------------------
    % Array shape: 'linear' or 'circular'
    ArrayShape = 'linear';
    switch ArrayShape
        case 'linear'
            arr_type = 1;
        case 'circular'
            arr_type = 2; %'circular' not supported in this code
    end
    %--------------------------------------------------------------------------
    % Number of microphones
    N_transducers = cfg.nmic;
    %--------------------------------------------------------------------------
    % Spacing of the microphones
    % - in case of linear array:      Spacing defines distance between microphones in meters
    % - in case of circular array:    Spacing defines radius of array in meters
    Spacing = 0; %if 0, arbitrary spacing is defined in BF_Array_Geometry.m            
    %--------------------------------------------------------------------------                                
    % constraint of the white noise gain (WNG) in dB
    % Note that the maximum possible value for the WNG is 10*log10(N_transducer),
    % which corresponds to delay-and-sum beamforming.
    Limit_WNG_dB =  wng_limit_db;%10*log10(N_transducers);
    %--------------------------------------------------------------------------  
    %Order of polynomial beamformer (P=0: normal FSB, P>0: polynomial
    %beamformer
    P=0;
    %--------------------------------------------------------------------------
    % desired look direction in degree (azimuth and elevation according to 
    % [Van Trees, Optimum Array Processing])
    LookDirection.azimuth = look_azimuth;
    LookDirection.elevation = look_elevation;%90 - atand(0.73/1.0); %(heigth of source w.r.t. array / distance in x-direction of source w.r.t. array)
    %--------------------------------------------------------------------------


    %--------------------------------------------------------------------------
    % Filter design -> call RobustFSB(...)
    %--------------------------------------------------------------------------
    [FilterImpulseResponse, PlantMatrix, BeamformerResponse, FrequencyAxis, AngularAxis, WhiteNoiseGain_dB] =...
            RobustFSBORG(N_transducers, Spacing, Limit_WNG_dB, P, 1, 1, arr_type, LookDirection, design);
    %--------------------------------------------------------------------------
end