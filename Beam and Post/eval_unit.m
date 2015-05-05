function [varargout] = eval_unit(speech_in, noise_in,N,varargin)
% [ssnre,lar,sd] = eval_unit(speech_in, noise_in,N,refsig,speech_out,noise_out,signal_out,porder)
% The Evaluation Unit calculates all possible objective measures
% (see diploma thesis, section 5.3)
% ssnre
% lar
% sd SSNR Enhancement SSNRE = SSNR_out - SSNR_in
% LAR (log-area-ratio distance) (optional)
% SD (speech degradation) (optional)
% speech_in
% noise_in
% N
% refsig
% speech_out
% noise_out
% signal_out
% porder input speech only (vector!)
% input noise only (vector!)
% Block length; default N = 256
% reference signal to calc. LAR and SD (vector!)
% output speech only (vector!)
% output noise only (vector!)
% output signal (vector!)
% model order for Levinson-Durbin Recursion (see thesis, appendix B.2);
% typ. porder = 14
% hint: to provide an accurate objective measures use signals that hardly contain
% intervals of silence in the speech utterances
% used functions: calc_parcor.m
% Syntax:
% to calculate SSNR of a signal
% [ssnr] = eval_unit(speech_in,noise_in,N)
% to calc. SSNRE, LAR & SD between the input signal and the processed signal
% [ssnre,lar,sd] = eval_unit(speech_in,noise_in,N,refsig,speech_out,noise_out,signal_out,porder)
if nargin<3 N=256; end
    if nargin<2
    help eval_unit
return;
end
% Declare input vectors
switch nargin
    case 3
    case 8
        refsig=varargin{1};
        speech_out=varargin{2};
        noise_out=varargin{3};
        signal_out=varargin{4};
        porder=varargin{5};
        refsig=refsig(:);speech_out=speech_out(:); % columnvectors!
        noise_out=noise_out(:); signal_out=signal_out(:); % columnvectors!
    otherwise
    error('Wrong number of input arguments')
end
speech_in=speech_in(:); noise_in=noise_in(:); % columnvectors!
% Calc. the vector length as a multiple of the block length
Nx_orig=length(speech_in);
M=floor(Nx_orig/N);
Nx=M*N;
% Calc. the objective measures for each frame
l=1;
for i=1:N:(Nx-N)
    k1=i:(i+N-1);
    % Calculate input-SSNR for each frame
    ssnr=10*log10(sum(speech_in(k1).^2)./sum(noise_in(k1).^2));
    ssnr_in(l)=ssnr;
    if nargin==8
        % Calculate output-SSNR for each frame
        ssnr_out(l)=10*log10(sum(speech_out(k1).^2)./sum(noise_out(k1).^2));
        % Calculate PARCOR coefficients and Area coefficients
        p_refsig(:,l)=calc_parcor(refsig(k1),porder);
        p_signal_out(:,l)=calc_parcor(signal_out(k1),porder);
        p_speech_out(:,l)=calc_parcor(speech_out(k1),porder);
        g_refsig(:,l)=(1+p_refsig(:,l))./(1-p_refsig(:,l));
        g_signal_out(:,l)=(1+p_signal_out(:,l))./(1-p_signal_out(:,l));
        g_speech_out(:,l)=(1+p_speech_out(:,l))./(1-p_speech_out(:,l));
        % Calculate LAR and SD for each frame
        lar(l)=sqrt((1/porder)*sum((abs(20*log10(g_refsig(:,l)./g_signal_out(:,l)))).^2));
        sd(l)=sqrt((1/porder)*sum((abs(20*log10(g_refsig(:,l)./g_speech_out(:,l)))).^2));
    end
    l=l+1;
end
switch nargin
    case 3
    % Calc. the averaged SSNR of a signal over all frames
    varargout{1}=sum(ssnr_in)/length(ssnr_in);
case 8
    % Calc. the averaged input and output SSNR over all frames
    msnr_out=sum(ssnr_out)/length(ssnr_out);
    msnr_in=sum(ssnr_in)/length(ssnr_in);
    % Calc. the averaged SSNRE over all frames
    varargout{1}=msnr_out-msnr_in;
    % Calc. the averaged LAR & SD over 95% of all frames
    lar_5=sort(lar);
    mlar_5=lar_5(1:ceil(length(lar_5)*0.95));
    sd_5=sort(sd);
    msd_5=sd_5(1:ceil(length(sd_5)*0.95));
    varargout{2}=sum(mlar_5)/length(mlar_5);
    varargout{3}=sum(msd_5)/length(msd_5);
    otherwise
    error('Wrong number of input arguments!!')
end