%% Homework Set 8
% Sam Nazari
% EECE 7312
% 24-November-2013


%% Simulation parameters
N   = 1000;               %Samples
M   = 5;                  %Array elements in ULA
C   = 10;
snir = 20*log10(1/C)
U   = zeros(M,N);
Rhat= zeros(M); 
rng(0,'twister');

%% Target
phi_0 = 0;              %Target location (deg)
th_0  = pi*sin(phi_0);  %Target elec. angle 

% Construct source signal
s0    = 1;

% Construct the steering vector
A0    = [1 1 1 1 1]';

% Construct the signal
S0 = s0*A0;

%% Inteference
phi_1 = 30*pi/180;              %Inteference location 
th_1  = pi*sin(phi_1);          %Inteference elec angle

for k = 1:N+1
% Construct int. signal
a     = 2*pi.*rand(1,1)-pi;     %random phase shift
s1    = C*exp(i*a)';

% Construct the steering vector
A1 = [exp(-1*i*th_1.*[0:M-1])]';

% Combine
S1 = s1*A1;

%% Noise
en = randn(M,1)*2*0.316;
z  = randn(M,1)*2*0.316;

%% Create U
U(:,k) = S0 + S1 +[en+1i*z];
%U = S0 + S1 +[en+1i*z];
end

%% Sample/Theoretical covariance matrix
Ruu = A0*A0'+abs(C)^2*A1*A1'+0.316*eye(5,5);  % Theoretical
for k = 1:N+1
    Rhat = Rhat+U*U';
    Rhat = Rhat/N;
end

% Now lets see how close the theoretical is to the numerical..
d = norm((Ruu-Rhat),'fro')
% This is not bad.

% Form the theoretical MVDR beamformer..
Wopt    = A0'*inv(Ruu);
W       = A0'*inv(Rhat);

%% Plot
phi1 = pi.*sin(-pi:pi/50:pi);
for k=1:length(phi1)
    A = [exp(-1*i*phi1(:,k).*[0:M-1])]';
    pl(k)= W*A;
    popt(k)=Wopt*A;
end
p=linspace(-pi,pi,length(phi1))*180/pi;
figure,
plot(p,20*log10(abs(pl)),'--ks','LineWidth',1,'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]),hold on,
plot(p,20*log10(abs(popt)),'r--o','LineWidth',1,'MarkerEdgeColor','b',...
    'MarkerFaceColor',[1,1,.25]),grid on
legend('Beamforming using SMI','Theoretical (numerical)')
xlabel('Degrees'),ylabel('Magnitude in dB')
title('MVDR Beamformer Amplitude Response')