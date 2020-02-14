%% ME367 Project 1
% Author: Scott Dolan
% Date  : 02-05-2020
% This code solves for forces involved with driveshaft
clc; clear all; close all;
%% Solve for input torque T from given input power(P,hp) and rotational velocity(rpm)
P = 0.75; % hp
w = 600.; % rpm
T = P*550*12/(w*2*pi/60); % lb.in
dshaft = 1.2 % in
%% Solve for force F_B at input gear
thetaB = 20.; % deg
dB = 8.; % in
FBt = 2*T/dB; % lb
FB = FBt/cosd(thetaB); % lb
FBr = FB*sind(thetaB); % lb
fprintf('The output force at gear A and its tangential and radial components:\n')
fprintf('Force FB = %.3f lb\n',FB);
fprintf('Force FBr = %.3f lb\n',FBr);
fprintf('Force FBt = %.3f lb\n\n',FBt);

%% Solve for force FA at output gear

thetaA = 20.; % degrees;
dA = 20.; % in
FAt = 2*T/dA; % lb
FA = FAt/cosd(thetaA); % lb
FAr = FA*sind(thetaA); % lb
fprintf('The input force at gear A and its tangential and radial components:\n')
fprintf('Force FA = %.3f lb\n',FA);
fprintf('Force FAr = %.3f lb\n',FAr);
fprintf('Force FAt = %.3f lb\n\n',FAt);

% FA can also be solved using the 3D moment sum taken at origin
% As you can see, the values match

Fat = FBt*dB/dA;       % lb
Fa = Fat/cosd(thetaA); % lb
Far = Fa*sind(thetaA); % lb
%% Solve for reactions at both bearings
% dOA, dAC, dCD are the distances between the origin\first bearing & gear
% A, Gear A & 2nd bearing at C, and 2nd bearing at C and gear B.

dOA = 16.; % in
dAC = 14.; % in
dCB = 9.; % in
dtot = dOA+dAC+dCB; % in

fprintf('Shaft length from bearing at point O to gear A, dOA = %.1f\n',dOA);
fprintf('Shaft length from gear A to bearing at point C, dAC = %.1f\n',dAC);
fprintf('Shaft length from bearing at point C to gear B, dCB = %.1f\n',dCB);
fprintf('Total shaft length, dtotal = %.1f\n\n',dtot);

RCz = (FBt*(dtot)-FAr*dOA)/(dOA+dAC); % lb
RCy = (FBr*dtot-FAt*dOA)/(dOA+dAC); % lb
fprintf('Bearing at point C reaction in y-direction, RCy = %.3f\n',RCy);
fprintf('Bearing at point C reaction in z-direction, RCz = %.2f\n\n',RCz);

FD = -2000; % lb 
ROx = -FD;
ROy = -FBr+FAt+RCy; % lb
ROz = FAr-FBt+RCz; % lb
fprintf('Bearing at point O reaction in x-direction, ROx = %.1f\n',ROx);
fprintf('Bearing at point O reaction in y-direction, ROy = %.3f\n',ROy);
fprintf('Bearing at point O reaction in z-direction, ROz = %.3f\n\n',ROz);

fsum_y = -ROy + FAt+RCy-FBr; % lb
fsum_z = ROz-FAr-RCz+FBt; %lb

%%
x = 0:0.01:(dtot+dtot*.01);    % "

% Singularity functions for x-z plane
qxz = @(x) ROz .* (x==0) - FAr .* (x==dOA) - RCz.*(x==(dOA+dAC)) + FBt .* (x==dtot);
Vxz = @(x) ROz .* (x>=0).*x.^0 - FAr.*(x>=dOA).*(x-dOA).^0 - RCz .* (x>=(dOA+dAC)).*(x-(dOA+dAC)).^0 ...
    + FBt.*(x>=dtot).*(x-dtot).^0;
My  = @(x) ROz .* (x>=0).*x - FAr.*(x>=dOA).*(x-dOA) - RCz .* (x>=(dOA+dAC)).*(x-(dOA+dAC)) ...
    + FBt.*(x>=dtot).*(x-dtot);
 

figure(1);
subplot(3,1,1);
plot(x,qxz(x)); grid on;
xlabel('Distance along the x axis (in)')
ylabel('Force in x-z plane');
title('Force in x-z plane');

subplot(3,1,2);
plot(x,Vxz(x));
grid on;
xlabel('Distance along the z axis (in)')
ylabel('Shear Force in x-z plane');
title('Shear Force in x-z plane');

subplot(3,1,3);
plot(x,My(x));
xlabel('Distance along the z axis (in)')
ylabel('Moment in x-z plane (lb in)')
grid on
title('Moment in the x-z plane');

% Singularity functions for y-z plane
qxy = @(x) -ROy.*(x==0) + FAt.*(x==dOA) + RCy.*(x==(dOA+dAC)) - FBr.*(x==dtot);
Vxy = @(x) -ROy.*(x>=0).*x.^0 + FAt.*(x>=dOA).*(x-dOA).^0 + RCy.*(x>=(dOA+dAC)).*(x-(dOA+dAC)).^0 ...
    -FBr.*(x>=dtot).*(x-dtot).^0;
Mz = @(x) -ROy.*(x>=0).*x + FAt.*(x>=dOA).*(x-dOA) + RCy.*(x>=(dOA+dAC)).*(x-(dOA+dAC)) ...
    -FBr.*(x>=dtot).*(x-dtot);


figure(2)

subplot(3,1,1);
plot(x,qxy(x)); grid on;
xlabel('Distance along the x axis (in)')
ylabel('Force in x-y plane');
title('Force in x-y plane');

subplot(3,1,2);
plot(x,Vxy(x));
grid on;
xlabel('Distance along the x axis (in)')
ylabel('Shear Force in x-y plane');
title('Shear Force in x-y plane');

subplot(3,1,3);
plot(x,Mz(x))
xlabel('Distance along the x axis (in)')
ylabel('Moment in x-y plane (lb in)')
grid on
title('Moment in the x-y plane');

figure(3);
subplot(2,1,1);
Tx = @(x) 0.*(x<=dOA)+T.*(x>=dOA)-T.*(x>=dtot);
plot(x,Tx(x));
title('Torsion Diagram along x [lb-in]');
xlabel('Distance along the x axis (in)')
ylabel('Torsion [lb-in]');

subplot(2,1,2);
A = @(x) FD.*(x>=0)-FD.*(x>=dtot);
plot(x,A(x));
title('Axial Force Diagram [lb]');
xlabel('Distance along the x axis (in)')
ylabel('Axial Force [lb]');

figure(4);
subplot(3,1,1);
plot(x,sqrt(qxy(x).^2+qxz(x).^2));
xlabel('Distance along the x axis (in)')
ylabel('Force');
title('Magnitude of Force');

subplot(3,1,2);
plot(x,sqrt(Vxy(x).^2+Vxz(x).^2));
grid on;
xlabel('Distance along the x axis (in)')
ylabel('Shear Force');
title('Magnitude of Shear Force');

subplot(3,1,3);
plot(x,sqrt(Mz(x).^2 + My(x).^2));
xlabel('Distance along the x axis (in)')
ylabel('Total Moment Magnitude (lb in)')
grid on
title('Combined Moment Magnitude');


%%  6 obtain the bending moment at sections of interest
V = @(x) (sqrt(Vxy(x).^2+Vxz(x).^2));
M = @(x) (sqrt(Mz(x).^2 + My(x).^2));

x = dOA;       % in
VA = V(x);     % lb 
MA = M(x);     % lb-in
fprintf(1,'At x = %g in VA = %g lb\n',x,VA);
fprintf(1,'At x = %g in MA = %g lb-in\n\n',x,MA);

x = (dOA+dAC); % in
VC = V(x);     % lb
MC = M(x);     % lb-in
fprintf(1,'At x = %g in VC = %g lb\n',x,VC);
fprintf(1,'At x = %g in MC = %g lb-in\n\n',x,MC);

x = dtot;      % in
VB = V(x);     % lb
MB = M(x);     % lb-in
fprintf(1,'At x = %g in VB = %g lb\n',x,VB);
fprintf(1,'At x = %g in MB = %g lb-in\n\n',x,MB);
%% Finding Maximum Stresses in Shaft

r = dshaft*.5; % in 
I = .25*pi*r^4;
J = pi*.5*r^4;
Q = .5*pi*r^2*4*r/(3*pi); % in^3
b = dshaft; % in

% Maximum stresses in shaft from combined shear and moments
fprintf('Maximum stresses in shaft from combined shear and moments:\n');
sbMc = MC*r/I;
tau_max = VC*Q/(I*dshaft);
fprintf('Maximum transverse shear stress in shaft [%.2f in<= x <=%.2f in], tau_max = %.2f psi\n',dOA,dOA+dAC,tau_max);
fprintf('Maximum bending stress at point C, x = %.2f in, sigma_Mc = %.2f psi\n\n',dOA+dAC,sbMc);

% Point C at x = dOA + dAC = 30 in
fprintf('The maximum stresses at point C:\n\n')

% Stresses in xy-plane
fprintf('Stresses in xy-plane:\n');

% Transverse Shear Stress in shaft, tau_Vxy
Vxyc = Vxy(dOA+dAC); % lb
tau_Vxyc = Vxyc*Q/(I*b); % psi
fprintf('Transverse Shear Stress in shaft xy plane, tau_Vxyc = %.2f psi\n',tau_Vxyc);

% Transverse Shear Stress in shaft, tau_Vxz
Vxzc = Vxz(dOA+dAC); % lb
tau_Vxzc = Vxzc*Q/(I*b); % psi
% Sigma Bending, sbMc
Myc = My(dOA+dAC); % lb-in
Mzc = Mz(dOA+dAC); % lb-in
sigma_Myc = Myc*r/I;
sigma_Mzc = Mzc*r/I;


% Stress from applied torque, T, and axial force FD
% Shear from torsion in shaft, tau_Tc
tau_Tc = T*r/J; % psi
fprintf('Shear from torsion in shaft, tau_Tc = %.2f psi\n',tau_Tc);
sigma_FD = FD/(pi*r^2);
fprintf('Axial stress from applied force FD, sigma_FD = %.2f psi\n\n',sigma_FD);
%% Finding principal stresses at point C
%tau_xy = tau_Tc-tau_Vc;
C2 = sbMc;

%% Angle of Twist, Deflection, and Slope
E = 29E3; % ksi, Modulus of Elasticity
thetaz = @(x) (-ROy.*(x>=0).*(x.^2)./2 + FAt.*(x>=dOA).*((x-dOA).^2)./2 + RCy.*(x>=(dOA+dAC)).*((x-(dOA+dAC)).^2)./2 ...
    -FBr.*(x>=dtot).*((x-dtot).^2)./2)./(E*I);
Yz = @(x) (-ROy.*(x>=0).*(x.^3)./6 + FAt.*(x>=dOA).*((x-dOA).^3)./6 + RCy.*(x>=(dOA+dAC)).*((x-(dOA+dAC)).^3)./6 ...
    -FBr.*(x>=dtot).*((x-dtot).^3)./6)./(E*I);
figure(5);
plot(x,Yz(x));