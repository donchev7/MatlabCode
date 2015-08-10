function CDR = CDR_prop2(Cxx,Cnn,TDOA,cfg)
Css = exp(1i*2*pi*cfg.frange'*TDOA); 
% Cnn = besselj(0,2*pi*cfg.frange'*micspacing/cfg.c); 
Css = bsxfun(@times, ones(size(Cxx)), Css); 
Cnn = bsxfun(@times, ones(size(Cxx)), Cnn); 
% limit the magnitude of Cxx to prevent numerical problems
magnitude_threshold = 1-1e-10;
critical = abs(Cxx)>magnitude_threshold;
Cxx(critical) = magnitude_threshold .* Cxx(critical) ./ abs(Cxx(critical));
% apply CDR estimator
CDR = 1./(-abs(Cnn-exp(1j*angle(Css)))./(Cnn.*cos(angle(Css))-1)).*abs((exp(-1j*angle(Css)).*Cnn - (exp(-1i*angle(Css)).*Cxx))./(real(exp(-1i*angle(Css)).*Cxx) - 1));
% Ensure we don't get any negative or complex results due to numerical effects
CDR = max(real(CDR),0);
end