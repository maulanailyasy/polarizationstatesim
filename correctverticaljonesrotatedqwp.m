% Set the initial parameters
Jvertilinear = [0 ; 1];  % Jones vector for vertical linear polarization
Im = 1j;  % Imaginary unit
theta_deg = linspace(0, 180, 9);  % Define the rotation angles (for plotting the polarization state)
theta = deg2rad(theta_deg);  % Convert angles to radians
qwp = [1 0 ; 0 -Im];  % Jones matrix for a QWP

% Define rotation matrix
rot = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];

% Rotate the QWP and determine the resultant Jones matrix
Jafterqwp = zeros(2, length(theta));
for i = 1:length(theta)
    %rotated_qwp = rot(-theta(i))*qwp_rotated*rot(theta(i));
    rotated_qwp = rot(theta(i))*qwp*rot(-theta(i));
    Jafterqwp(:,i) = rotated_qwp*Jvertilinear;
end

% Extract magnitude and phase from Jones vectors
magx = abs(Jafterqwp(1,:));
phasex = angle(Jafterqwp(1,:));
magy = abs(Jafterqwp(2,:));
phasey = angle(Jafterqwp(2,:));

% Calculate the phase difference delta for each theta
delta = phasey - phasex;

% Calculate the components Ex and Ey for each theta and t
t = linspace(0, 2*pi, 1000);
Ex = zeros(length(t), length(theta));
Ey = zeros(length(t), length(theta));
for i = 1:length(theta)
    Ex(:,i) = magx(i)*cos(t);
    Ey(:,i) = magy(i)*(cos(delta(i))*cos(t) - sin(delta(i))*sin(t));
end

% Define the extracted colors, for color customization
colors = [
    0.6784, 0.8627, 0.1882;
    0.3686, 0.7882, 0.3843;
    0.1529, 0.6784, 0.5059;
    0.1294, 0.5608, 0.5529;
    0.1765, 0.4392, 0.5569;
    0.2392, 0.3059, 0.5412;
    0.2824, 0.1569, 0.4706;
    0.2824, 0.1373, 0.4549;
    103/255, 120/255, 182/255; % Light Blue
    51/255, 72/255, 152/255;  % Medium Blue
    153/255, 164/255, 206/255; % Sky Blue
    27/255, 41/255, 94/255;   % Dark Blue
    17/255, 7/255, 19/255;    % Almost Black
    237/255, 28/255, 7/255;   % Red
    126/255, 0/255, 4/255;    % Dark Red
    253/255, 181/255, 2/255;  % Orange
];

% Plot the polarization ellipses directly from Ex and Ey, with adjustments
figure;

for i = 1:length(theta)
    subplot(3, 3, i);  % Adjust the grid size as needed
    
    % Plot the ellipse
    plot(Ex(:,i), Ey(:,i), 'b', 'LineWidth', 1.0);
    grid on
    axis equal;
    xlim([-1.0 1.0]);
    ylim([-1.0 1.0]);
    %if abs(mag_comp1(i)) < 0.1
       %xlabel('E_x');
       %ylabel('E_y');
    %end
    
    
    title(['\theta = ' num2str(theta_deg(i)) '°']);
    
    % Add annotations for the polarization state
    if abs(delta(i) - pi/2) < 0.1 && abs(magx(i) - magy(i)) < 0.1
        text(0.6, 0.9, 'RC', 'Color', colors(6, :));
    elseif abs(delta(i) + pi/2) < 0.1 && abs(magx(i) - magy(i)) < 0.1
        text(0.6, 0.9, 'LC', 'Color', colors(8, :));
    elseif delta(i) > 0 && abs(magx(i)) > 0.1
        text(0.6, 0.9, 'RE', 'Color', colors(4, :));
    elseif delta(i) < 0 && abs(magx(i)) > 0.1
        text(0.6, 0.9, 'LE', 'Color', colors(4, :));
    elseif abs(magx(i)) < 0.1
        text(0.6, 0.9, 'V', 'Color', colors(2, :));
    end
end
%}


figure;

% Plotting the magnitude of both components
% Font size settings
titleFontSize = 11; % Adjust as needed
labelFontSize = 11;  % Adjust as needed

subplot(2,2,1);
plot(theta_deg, magx, 'b', 'LineWidth', 1.0);
title('J_x component ', 'FontSize', titleFontSize);
%xlabel('\theta (deg)', 'FontSize', labelFontSize);
ylabel('Magnitude', 'FontSize', labelFontSize);

subplot(2,2,2);
plot(theta_deg, magy, 'r', 'LineWidth', 1.0);
title('J_y component');
%xlabel('\theta (deg)');
ylabel('Magnitude');

% Plotting the phase of both components
subplot(2,2,3);
plot(theta_deg, rad2deg(phasex), 'b', 'LineWidth', 1.0);
%title('Phase of Component J_x vs. \theta (deg)');
xlabel('\theta (deg)');
ylabel('Phase (degrees)');
%xlim([-1.0 1.0]);

subplot(2,2,4);
plot(theta_deg, rad2deg(phasey), 'r', 'LineWidth', 1.0);
%title('Phase of Component J_y vs. \theta (deg)');
xlabel('\theta (deg)');
ylabel('Phase (degrees)');