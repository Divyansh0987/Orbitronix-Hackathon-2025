
% GOAL:
%   Model and visualize a satellite's orbit from given orbital parameters 
%   (RAAN, inclination, radius, orbital period, argument of perigee) 
%   and plot in spherical coordinates.
%
% STEPS:
%   1) Define orbital elements for a circular orbit
%   2) Generate perifocal coordinates
%   3) Rotate to ECI
%   4) Convert ECI to spherical (r, theta, phi)
%   5) Plot spherical coordinates vs. time and also a 3D orbit

clear; clc; close all;

% 1) DEFINE ORBITAL PARAMETERS (CIRCULAR ORBIT)
mu_Earth = 3.986004418e5;  % [km^3/s^2] Earth standard gravitational parameter
R_Earth  = 6371;           % [km] approximate Earth radius for reference

%define time in seconds
hour = 3;
totalTime = hour*3600;

% customised orbit elements:
incl_deg = 10;  % Inclination i [deg]
raan_deg = 72;  % Right Ascension of Ascending Node Ω [deg]
w_deg = 120.479;  % Argument of Perigee ω [deg]
r_orbit = 6771.6;  % [km] for a circular orbit (constant radius)
T_orbit = 2*pi * sqrt(r_orbit^3 / mu_Earth);  
n = 2*pi / T_orbit;  % Mean motion [rad/s]

% aqua parameterss
aqua_i = deg2rad(98.199);  % Inclination i [deg]
aqua_raan = deg2rad(95.206);  % Right Ascension of Ascending Node Ω [deg]
aqua_w = deg2rad(120.479);  % Argument of Perigee ω [deg]
aqua_r = 7080.6;  % [km] for a circular orbit (constant radius)
aqua_T = 2*pi * sqrt(aqua_r^3 / mu_Earth);  
aqua_n = 2*pi / aqua_T;  % Mean motion [rad/s]

% convert angles to radians
i_rad       = deg2rad(incl_deg);
raan_rad   = deg2rad(raan_deg);
w_rad       = deg2rad(w_deg);

% some orbital info
disp('Aqua Orbit radius: 7080.6 km')
disp('Orbital period: 5929.5 s (1.6 hours)')
fprintf('Orbit Radius: %.1f km\n', r_orbit);
fprintf('Orbital period: %.1f s (%.1f hours)\n', T_orbit, T_orbit/3600);

% 2) TRUE ANOMALY & PERIFOCAL COORDINATES
% parameterize the orbit by true anomaly nu in [0, 2π].
% For a circular orbit (eccentricity e=0), the radius is constant
numPoints = totalTime; % number of steps in one orbit
t_vec = linspace(0, totalTime, numPoints);  % optional time vector
nu_vec = n * t_vec;% if starting at nu=0 at t=0
nu_vec_aqua = aqua_n * t_vec;

% ensure nu_vec goes from 0 to 2π exactly:
nu_vec = mod(nu_vec, 2*pi); 
nu_vec_aqua = mod(nu_vec_aqua, 2*pi); 

% perifocal frame position for a circular orbit:
r_pf = [ r_orbit * cos(nu_vec);
         r_orbit * sin(nu_vec);
         zeros(1,numPoints) ];
r_pf_a = [ aqua_r * cos(nu_vec_aqua);
         aqua_r * sin(nu_vec_aqua);
         zeros(1,numPoints) ];

% 3) 
% Standard rotation sequence from perifocal to ECI:
% r_ECI = Rz(-Ω) * Rx(-i) * Rz(-ω) * r_pf

% Define rotation matrices
Rz_minusOmega = [ cos(-raan_rad),  sin(-raan_rad), 0
                 -sin(-raan_rad),  cos(-raan_rad), 0
                  0,                0,               1 ];
Rx_minusi     = [ 1,  0,           0
                  0,  cos(-i_rad),  sin(-i_rad)
                  0, -sin(-i_rad),  cos(-i_rad) ];
Rz_minusw     = [ cos(-w_rad),  sin(-w_rad), 0
                 -sin(-w_rad),  cos(-w_rad), 0
                  0,            0,           1 ];
Rz_minusOmega_a = [ cos(-aqua_raan),  sin(-aqua_raan), 0
                 -sin(-aqua_raan),  cos(-aqua_raan), 0
                  0,                0,               1 ];
Rx_minusi_a     = [ 1,  0,           0
                  0,  cos(-aqua_i),  sin(-aqua_i)
                  0, -sin(-aqua_i),  cos(-aqua_i) ];
Rz_minusw_a     = [ cos(-aqua_w),  sin(-aqua_w), 0
                 -sin(-aqua_w),  cos(-aqua_w), 0
                  0,            0,           1 ];
% Combined rotation matrix
Q = Rz_minusOmega * Rx_minusi * Rz_minusw;
Qa = Rz_minusOmega_a * Rx_minusi_a * Rz_minusw_a;


% transform all points
r_eci = Q * r_pf;
r_eci_a = Qa * r_pf_a;
x_eci = r_eci(1,:);
y_eci = r_eci(2,:);
z_eci = r_eci(3,:);
x_eci_a = r_eci_a(1,:);
y_eci_a = r_eci_a(2,:);
z_eci_a = r_eci_a(3,:);

% 4) convert ECI to spherical coordinates
r_sph     = sqrt(x_eci.^2 + y_eci.^2 + z_eci.^2);
theta_sph = atan2( sqrt(x_eci.^2 + y_eci.^2), z_eci ); % -> polar angle from +Z axis
phi_sph   = atan2( y_eci, x_eci ); % -> azimuth angle in XY plane

r_sph_a     = sqrt(x_eci_a.^2 + y_eci_a.^2 + z_eci_a.^2);
theta_sph_a = atan2( sqrt(x_eci_a.^2 + y_eci_a.^2), z_eci_a );
phi_sph_a   = atan2( y_eci_a, x_eci_a ); 


% (A) plot spherical coords vs. time (or vs. anomaly)
figure('Name','Spherical Coordinates','Color','w');

subplot(2,1,1);
plot(linspace(0, 18000, numPoints), rad2deg(theta_sph), 'r', 'LineWidth',1.5);
hold on;
plot(linspace(0, 18000, numPoints), rad2deg(theta_sph_a), 'b', 'LineWidth',1.5);
xlabel('Time [s]');
ylabel('\theta [deg]');
title('Polar Angle vs. Time');
grid on;

subplot(2,1,2);
plot(linspace(0, 1800, numPoints), rad2deg(phi_sph), 'r', 'LineWidth',1.5);
hold on;
plot(linspace(0, 1800, numPoints), rad2deg(phi_sph_a), 'b', 'LineWidth',1.5);
xlabel('Time [s]');
ylabel('\phi [deg]');
title('Azimuth Angle vs. Time');
grid on;

% (B) 3D Orbit Plot in ECI (for reference)
figure('Name','3D Orbit in ECI','Color','w');
plot3(x_eci, y_eci, z_eci, 'r', 'LineWidth',1.5);
hold on; grid on; axis equal;
plot3(x_eci_a, y_eci_a, z_eci_a, 'b', 'LineWidth',1.5);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Satellite Orbit in ECI Coordinates');

% Draw Earth for reference
[XX, YY, ZZ] = sphere(30);
surf(R_Earth*XX, R_Earth*YY, R_Earth*ZZ, ...
    'FaceColor',[0.5 0.7 1], 'EdgeColor','none','FaceAlpha',0.5);

legend('Custom Orbit','Aqua Orbit','Location','best');