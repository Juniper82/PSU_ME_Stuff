%% ICE 1
% Shaft Design for Steady Torsion and Fully Reversed Bending
% Revision: Spring 2019
 clear; clc; close all;
%sets the default line width in plots.  Needs to be done
%only once, but it is done here so we are sure.
set(0,'DefaultLineLineWidth',2);

%% some conversions defined
hp_to_lb_in_per_sec = 6600; % 1 hp = 6600 lb*in/sec
rpm_to_radper_sec = 2*pi/60; % 1 rmp = 2*pi/60 rads/sec
deg_to_rad = (2*pi / 360); % conversion from deg to radians

%% Given
P = 2;          % hp
OMEGA = 1725;   % rpm
% Starting guessed values for some of the shaft dimensions
a = 1.5;        % in
b = 5.0;        % in
c = 6.5;        % in
p = 2.0;        % in
q = 6.75;       % in
s = 1.50;       % in
dg = 6.00;      % in
rg = 0.5*dg;     % in
ds = 6.00;      % in
rs = 0.5*ds;     % in
fprintf(1,'Gear Radius = %g in, Sheave Radisu = %g in \n',rg,rs);
% Gear Pressure Angle
PHIdeg = 20;    % deg
PHI = PHIdeg * deg_to_rad;   % rad
fprintf(1,'PHI = %g radians\n',PHI);
%% Solution
%% 1 Find the toruqe based on input power
T = P*hp_to_lb_in_per_sec/(OMEGA * rpm_to_radper_sec);     % lb in
fprintf(1,'T = %g lb*in\n',T);

%% 2 Find the belt force and total belt force
% Using the above Torque and belt tension ratio as 5
RHO = 5;     % Tension ratio - no units
Fn = T/rs;   % lb
F1 = Fn / (1 - (1/RHO)); % lb
F2 = F1 - Fn;    % lb
Fs = F1 + F2 ;   % lb
fprintf(1,'Fn = %g lb, F1 = %g lb\n',Fn,F1);
fprintf(1,'F2 = %g lb, Fs = %g lb\n',F2,Fs);
%
%% 3 Find the force due to the gear
Fgt = T / rg;    % lb
Fgr = Fgt * tan(PHI);    % lb
fprintf(1,'Fgt = %g lb, Fgr = %g lb\n',Fgt,Fgr);
%
%% 4 Assuming gear and sheave force to be concentrated at their center
% remember that these are distributed forces in reality
% Below reactions are based on the FBD drawn
% 
R2x = (1/b)*(-Fgr*p - Fs*q);  % lb
R1x = -Fgr - R2x - Fs;   % lb
R2y = (1/b)*(Fgt*p); % lb
R1y = Fgt - R2y;     % lb
fprintf(1,'R1x = %g lb, R1y = %g lb, \nR2x = %g lb, R2y = %g lb\n',R1x,R1y,R2x,R2y);
%
%% 5 Discontunity functions and obtaining the shear force and bending moments
% using the discontinuity function (as shown in the handout)
% we will then plot the loading function, shear force and bending moment in
% each place and then finally the total magnitude

% define the load function
z = 0:0.01:7;    % in
qxz = @(z) R1x .* (z==0) + Fgr .* (z==p) + R2x .*(z==b) + Fs .* (z==q);
Vxz = @(z) R1x .* (z>=0).* z.^0 + Fgr .* (z>=p).*(z-p).^0 + ...
    R2x .*(z>=b).* (z-b).^0 + Fs .* (z>=q).*(z-q).^0;
My = @(z) R1x .* (z>=0) .* z + Fgr .* (z>=p) .*...
    (z-p) + R2x .* (z>=b) .* (z-b) + Fs .* (z>=q) .* (z-q);

figure(1);
subplot(3,1,1);
plot(z,qxz(z)); grid on;
xlabel('Distance along the z axis (in)');
ylabel('Force in x-z plane');
title('Force in x-z plane');

subplot(3,1,2);
plot(z,Vxz(z));
grid on;
xlabel('Distance along the z axis (in)');
ylabel('Shear Force in x-z plane');
title('Shear Force in x-z plane');

subplot(3,1,3);
plot(z,My(z));
xlabel('Distance along the z axis (in)');
ylabel('Moment in x-z plane (lb in)');
grid on
title('Moment in the x-z plane');

% similarily compute the forces in the yz place

qyz = @(z) R1y .* (z==0) - Fgt .* (z==p) + R2y .* (z==b);
Vyz = @(z) R1y .*(z>=0).* z.^0 - Fgt .*(z>=p).* (z-p).^0 + ...
    R2y .*(z>=b).* (z-b).^0;
Mx = @(z) R1y .* (z>=0) .* z - Fgt .* (z>=p) .*...
     (z-p) + R2y .* (z>=b) .* (z-b);

figure(2)

subplot(3,1,1);
plot(z,qyz(z)); grid on;
xlabel('Distance along the z axis (in)');
ylabel('Force in y-z plane');
title('Force in y-z plane');

subplot(3,1,2);
plot(z,Vyz(z));
grid on;
xlabel('Distance along the z axis (in)');
ylabel('Shear Force in y-z plane');
title('Shear Force in y-z plane');

subplot(3,1,3);
plot(z,Mx(z))
xlabel('Distance along the z axis (in)');
ylabel('Moment in y-z plane (lb in)');
grid on
title('Moment in the y-z plane');

figure(3)

subplot(3,1,1); 
% magnitude of the load along the shaft
Q = @(z) ((qxz(z).^2 + qyz(z).^2).^(1/2));
plot(z,Q(z));
grid on;
xlabel('Distance along the z axis (in)');
ylabel('Force');
title('Magnitude of Force');

% magnitude of the shear force along the shaft
subplot(3,1,2);
V = @(z)(Vxz(z).^2 + Vyz(z).^2).^(1/2);
plot(z,V(z));
grid on;
xlabel('Distance along the z axis (in)');
ylabel('Shear Force');
title('Magnitude of Shear Force');

% magnitude of the bending moment along the shaft
subplot(3,1,3);
M = @(z)(Mx(z).^2 + My(z).^2).^(1/2);
plot(z,M(z))
xlabel('Distance along the z axis (in)')
ylabel('Total Moment magnitude (lb in)')
grid on
title('Combined Moment magnitude');
%%  6 obtain the bending moment at sections of interest
z = p   ;    % in
MB = M(z);   % lb in
fprintf(1,'At z = %g in MB = %g lb *in\n',z,MB);
z = b  ;     % in
MC = M(z);   % lb in
fprintf(1,'At z = %g in MC = %g lb *in\n',z,MC);
z = c ;      % in
MD = M(z);   % lb in
fprintf(1,'At z = %g in MD = %g lb *in\n',z,MD);

%Show points of interest on graph
Bhandle = text(p,MB,'B');
Chandle = text(b,MC,'C');
Dhandle = text(c,MD,'D');
