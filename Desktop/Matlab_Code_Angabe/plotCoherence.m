function Cxx = plotCoherence(noise,alpha,mic_distance,frange,theoretical_coherence)
%theoretical_coherence
% Measurement of the coherence function between two channels of a recorded
% multi-channel noise (late tail in RIR modeled as noise)
% and comparison with the ideal sin(x)/x - coherence function
Pxx = estimate_psd(noise,alpha);
nmic = length(mic_distance);
Cxx=zeros(length(frange),1);
Gamma=zeros(length(frange),1);
beta = (2*pi*frange)/342;
for m = 1:(nmic-1)
    for n = (m+1):nmic
        Cxx = Cxx + estimate_cpsd(noise(:,m),noise(:,n),alpha)./sqrt(Pxx(:,m).*Pxx(:,n));
        if strcmp(theoretical_coherence,'sinc')
            Gamma = Gamma + sinc(beta*mic_distance(m,n));
        elseif strcmp(theoretical_coherence,'bessel')
            Gamma = Gamma + besselj(0,beta*mic_distance(m,n));
        end
    end
end

Combination = 2/(nmic*(nmic-1));

Cxx=real(Cxx)*Combination;
Gamma=real(Gamma)*Combination;

plot(frange,real(Cxx),'r');
hold on
plot(frange,real(Gamma));
hold off
axis([frange(1) frange(end) -0.5 1])
legend('Measured',['Theory'],3)
xlabel('Frequency [Hz]')
ylabel(['Real(\Gamma)'])

end