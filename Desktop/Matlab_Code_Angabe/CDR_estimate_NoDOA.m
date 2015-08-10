function W = CDR_estimate_NoDOA(Cxx,Cnn)
Cnn = bsxfun(@times, ones(size(Cxx)), Cnn); % extend to dimension of Cxx
% limit the magnitude of Cxx to prevent numerical problems
magnitude_threshold = 1-1e-10;
critical = abs(Cxx)>magnitude_threshold;
Cxx(critical) = magnitude_threshold .* Cxx(critical) ./ abs(Cxx(critical));

W =  (-(abs(Cxx).^2 + Cnn.^2.*real(Cxx).^2 - Cnn.^2.*abs(Cxx).^2 - 2.*Cnn.*real(Cxx) + Cnn.^2).^(1/2) - abs(Cxx).^2 + Cnn.*real(Cxx))./(abs(Cxx).^2-1);

% Ensure we don't get any negative or complex results due to numerical effects
W = max(real(W),0);
end
