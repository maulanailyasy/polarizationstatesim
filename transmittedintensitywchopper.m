% Initial Jones vector for the polarization state of light/laser beam
VP = [0; 1];
lp = [1 0; 0 0];
retro = [-1 0; 0 -1];
M1 = [1 0; 0 1];
M2 = [1 0; 0 -1]; % Assume M2 = M3 = M4 are the same
SM1 = [-1 0; 0 -1];
alphadeg = 0;
alpha = deg2rad(alphadeg);
rot = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
rotSM = rot*SM1;

% Define analyzer's Jones matrix
phi = 90; % Angle of analyzer
phirad = deg2rad(phi);
rotaposLP = [cos(phirad) -sin(phirad); sin(phirad) cos(phirad)];
rotanegLP = [cos(phirad) -sin(-phirad); sin(-phirad) cos(phirad)];
LP = rotaposLP*lp*rotanegLP;

% Define new vectors for time
%times = linspace(0, 0.06, 17500000); 
times = linspace(0, 0.06, 1000); %for faster simulation
intensity_simulated = zeros(1, length(times));

% Define quarter-wave plate's Jones matrix with synchronized phase
theta_base_freq = 57; % frequency of quarter-wave plate rotation in Hz
theta_phase = @(time) 360*theta_base_freq*time; % function of time to theta phase in degrees

% Define quarter-wave plate's Jones matrix
Q = @(theta) [cosd(theta)^2 - 1i*sind(theta)^2, -(1 + 1i)*cosd(theta)*sind(theta);
             -(1 + 1i)*cosd(theta)*sind(theta), sind(theta)^2 - 1i*cosd(theta)^2];

% Chopper function for vertical polarization
chop_freq = 114; % frequency of chopper in Hz
chop_period = 1/chop_freq; % period of chopper in seconds
chopper = @(time) double(mod(time, 1 / chop_freq) < (1 / chop_freq) / 2);

A = retro*VP;
B = M1*A;
for idx = 1:length(times)
    time = times(idx); % synchronized time for each sample
    theta = theta_phase(time); % synchronized phase of quarter-wave plate
    chopper_state = chopper(time); % get chopper state at this time
    C = Q(theta)*B;
    D = M2*C;
    E = M2*D;
    F = SM1*E;
    G = M2*F;
    H = rotSM*G;
    J = LP*H;
    % Multiply intensity by chopper state (0 or 1)
    intensity_simulated(idx) = chopper_state*(abs(conj(J')*J));
end

figure;
plot(times, intensity_simulated, 'LineWidth', 1); % Plot simulation
%hold on;
%plot(times, normintensityofphase196deg, 'LineWidth', 1); % the experimental data
xlabel('Time (ms)');
ylabel('Transmitted Intensity');
legend('Simulated', 'Measured');
grid on;