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

LOA = 16.; % in
LAC = 14.; % in
LCB = 9.; % in
Ltot = LOA+LAC+LCB; % in

fprintf('Shaft length from bearing at point O to gear A, dOA = %.1f\n',LOA);
fprintf('Shaft length from gear A to bearing at point C, dAC = %.1f\n',LAC);
fprintf('Shaft length from bearing at point C to gear B, dCB = %.1f\n',LCB);
fprintf('Total shaft length, dtotal = %.1f\n\n',Ltot);

RCz = (FBt*(Ltot)-FAr*LOA)/(LOA+LAC); % lb
RCy = (FBr*Ltot-FAt*LOA)/(LOA+LAC); % lb
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
x = 0:0.01:(Ltot+Ltot*.01);    % "

% Singularity functions for x-z plane
qxz = @(x) ROz .* (x==0) - FAr .* (x==LOA) - RCz.*(x==(LOA+LAC)) + FBt .* (x==Ltot);
Vxz = @(x) ROz .* (x>=0).*x.^0 - FAr.*(x>=LOA).*(x-LOA).^0 - RCz .* (x>=(LOA+LAC)).*(x-(LOA+LAC)).^0 ...
    + FBt.*(x>=Ltot).*(x-Ltot).^0;
My  = @(x) ROz .* (x>=0).*x - FAr.*(x>=LOA).*(x-LOA) - RCz .* (x>=(LOA+LAC)).*(x-(LOA+LAC)) ...
    + FBt.*(x>=Ltot).*(x-Ltot);
 

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

% Singularity functions for x-y plane
qxy = @(x) -ROy.*(x==0) + FAt.*(x==LOA) + RCy.*(x==(LOA+LAC)) - FBr.*(x==Ltot);
Vxy = @(x) -ROy.*(x>=0).*x.^0 + FAt.*(x>=LOA).*(x-LOA).^0 + RCy.*(x>=(LOA+LAC)).*(x-(LOA+LAC)).^0 ...
    -FBr.*(x>=Ltot).*(x-Ltot).^0;
Mz = @(x) -ROy.*(x>=0).*x + FAt.*(x>=LOA).*(x-LOA) + RCy.*(x>=(LOA+LAC)).*(x-(LOA+LAC)) ...
    -FBr.*(x>=Ltot).*(x-Ltot);


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
Tx = @(x) 0.*(x<=LOA)+T.*(x>=LOA)-T.*(x>=Ltot);
plot(x,Tx(x));
title('Torsion Diagram along x [lb-in]');
xlabel('Distance along the x axis (in)')
ylabel('Torsion [lb-in]');

subplot(2,1,2);
A = @(x) FD.*(x>=0)-FD.*(x>=Ltot);
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

x = LOA;       % in
VA = V(x);     % lb 
MA = M(x);     % lb-in
fprintf(1,'At x = %g in VA = %g lb\n',x,VA);
fprintf(1,'At x = %g in MA = %g lb-in\n\n',x,MA);

x = (LOA+LAC); % in
VC = V(x);     % lb
MC = M(x);     % lb-in
fprintf(1,'At x = %g in VC = %g lb\n',x,VC);
fprintf(1,'At x = %g in MC = %g lb-in\n\n',x,MC);

x = Ltot;      % in
VB = V(x);     % lb
MB = M(x);     % lb-in
fprintf(1,'At x = %g in VB = %g lb\n',x,VB);
fprintf(1,'At x = %g in MB = %g lb-in\n\n',x,MB);
%% Finding Maximum Stresses in Shaft
% since diameter of shaft is constant I & J will be constant and for now we
% will not write I and J as a function of d
r = dshaft*.5; % in 
I = .25*pi*r^4;
J = pi*.5*r^4;
Q = .5*pi*r^2*4*r/(3*pi); % in^3
b = dshaft; % in

% Maximum stresses in shaft from combined shear and moments
fprintf('Maximum stresses in shaft from combined shear and moments:\n');
sbMc = MC*r/I;
tau_max = VC*Q/(I*dshaft);
fprintf('Maximum transverse shear stress in shaft [%.2f in<= x <=%.2f in], tau_max = %.2f psi\n',LOA,LOA+LAC,tau_max);
fprintf('Maximum bending stress at point C, x = %.2f in, sigma_Mc = %.2f psi\n\n',LOA+LAC,sbMc);

% Point C at x = dOA + dAC = 30 in
fprintf('The maximum stresses at point C:\n\n')

% Stresses in xy-plane
fprintf('Stresses in xy-plane:\n');

% Transverse Shear Stress in shaft, tau_Vxy
Vxyc = Vxy(LOA+LAC); % lb
tau_Vxyc = Vxyc*Q/(I*b); % psi
fprintf('Transverse Shear Stress in shaft xy plane, tau_Vxyc = %.2f psi\n',tau_Vxyc);

% Transverse Shear Stress in shaft, tau_Vxz
Vxzc = Vxz(LOA+LAC); % lb
tau_Vxzc = Vxzc*Q/(I*b); % psi
% Sigma Bending, sbMc
Myc = My(LOA+LAC); % lb-in
Mzc = Mz(LOA+LAC); % lb-in
sigma_Myc = Myc*r/I;
sigma_Mzc = Mzc*r/I;


% Stress from applied torque, T, and axial force FD
% Shear from torsion in shaft, tau_Tc
tau_Txc = T*r/J; % psi
fprintf('Shear from torsion in shaft, tau_Tc = %.2f psi\n',tau_Txc);
sigma_FD = FD/(pi*r^2);
fprintf('Axial stress from applied force FD, sigma_FD = %.2f psi\n\n',sigma_FD);
%% Finding principal stresses at point C on shaft
c2 = @(sigX,sigY,sigZ) sigX + sigY + sigZ;
c1 = @(tauXY,tauXZ,tauYZ,sigX,sigY,sigZ) tauXY^2 + tauYZ^2 + tauXZ^2 -sigX*sigY - sigY*sigZ - sigZ*sigX;
c0 = @(tauXY,tauXZ,tauYZ,sigX,sigY,sigZ) sigX*sigY*sigZ+2*tauXY*tauYZ*tauXZ-sigX*tauXZ^2-sigZ*tauXY^2;
 
% point A on cross section:
sigma_xA = sigma_FD-sigma_Myc;
tauYZa = 0;
tauXZa = 0;
tauXYa = tau_Txc + tau_Vxyc;
sigXa = sigma_xA;
sigYa = 0;
sigZa = 0;
c2a = c2(sigXa,sigYa,sigZa);
c1a = c1(tauXYa,tauXZa,tauYZa,sigXa,sigYa,sigZa);
c0a = c0(tauXYa,tauXZa,tauYZa,sigXa,sigYa,sigZa);
principalA = roots([1,c2a,c1a,c1a]);
sigma2 = getsigma2(principalA);
fprintf('The principal stresses at point A on cross section are:\n');
fprintf('sigma1 = %.2f\n',max(principalA));
fprintf('sigma2 = %.2f\n',sigma2);
fprintf('sigma3 = %.2f\n\n',min(principalA));

% point B on cross section:
sigma_xB = sigma_FD-sigma_Mzc;
tauXYb = 0;
tauYZb = 0;
tauXZb = tau_Txc + tau_Vxzc;
sigXb = sigma_xB;
sigYb = 0;
sigZb = 0;
c2b = c2(sigXb,sigYb,sigZb);
c1b = c1(tauXYb,tauXZb,tauYZb,sigXb,sigYb,sigZb);
c0b = c0(tauXYb,tauXZb,tauYZb,sigXb,sigYb,sigZb);
principalB = roots([1,c2b,c1b,c1b]);
sigma2 = getsigma2(principalB);
fprintf('The principal stresses at point B on cross section are:\n');
fprintf('sigma1 = %.2f\n',max(principalB));
fprintf('sigma2 = %.2f\n',sigma2);
fprintf('sigma3 = %.2f\n\n',min(principalB));

% point C on cross section:
sigma_xC = sigma_FD+sigma_Myc;
tauXYc = tau_Txc - tau_Vxyc;
tauYZc = 0;
tauXZc = 0;
sigXc = sigma_xC;
sigYc = 0;
sigZc = 0;
c2c = c2(sigXc,sigYc,sigZc);
c1c = c1(tauXYc,tauXZc,tauYZc,sigXc,sigYc,sigZc);
c0c = c0(tauXYc,tauXZc,tauYZc,sigXc,sigYc,sigZc);
principalC = roots([1,c2c,c1c,c1c]);
sigma2 = getsigma2(principalC);
fprintf('The principal stresses at point C on cross section are:\n');
fprintf('sigma1 = %.2f\n',max(principalC));
fprintf('sigma2 = %.2f\n',sigma2);
fprintf('sigma3 = %.2f\n\n',min(principalC));

% point D on cross section:
sigma_xD = sigma_FD+sigma_Mzc;
tauXYd = 0;
tauYZd = 0;
tauXZd = tau_Txc + tau_Vxzc;
sigXd = sigma_xD;
sigYd = 0;
sigZd = 0;
c2d = c2(sigXd,sigYd,sigZd);
c1d = c1(tauXYd,tauXZd,tauYZd,sigXd,sigYd,sigZd);
c0d = c0(tauXYd,tauXZd,tauYZd,sigXd,sigYd,sigZd);
principalD = roots([1,c2d,c1d,c1d]);
sigma2 = getsigma2(principalD);
fprintf('The principal stresses at point D on cross section are:\n');
fprintf('sigma1 = %.2f\n',max(principalD));
fprintf('sigma2 = %.2f\n',sigma2);
fprintf('sigma3 = %.2f\n\n',min(principalD));

%% Angle of Twist, Deflection, and Slope
% since diameter of shaft is constant I & J will be constant and for now we
% will not write I and J as a function of d
% Angle of twist
d = 1.2;
J = @(d) (pi/32)*d.^4; % area polar moment of inertia in inches ^4
I = @(d) (pi/64)*d.^4; % second area moment of inertia in inches ^4
%I = .25*pi*r^4;
%J = pi*.5*r^4;
x = 0:0.01:(Ltot+Ltot*.01);    % "
G = 11600E3; %ksi
E = 29000E3; % ksi, Modulus of Elasticity
Loa = LOA;
Lac = LAC;
Lcb = LCB;
Toa = 0;
Tac = T;
Tcb = T;
Tab = T;
Phi_oa = (Toa*Loa/(G*J(d)))*180/pi;
Phi_ac = (Tac*Lac/(G*J(d)))*180/pi;
Phi_cb = (Tcb*Lcb/(G*J(d)))*180/pi;
Phi = @(x) sum((Toa.*x/(G*J(d)).*(x==0) + Tab.*(x-LOA)./(G*J(d)).*(x>LOA)).*(180/pi));
Phi_shaft = Phi_oa + Phi_ac + Phi_cb;
fprintf('The total twist in the shaft, Phi_shaft = %g deg\n\n',Phi_shaft);
Phi_ab = Tab*(Lac+Lcb)/(G*J(d))*180/pi;
C4 = 0;
C3 = (ROy*(LOA+LAC)^3-FAt*LAC^3)/(6*(LOA+LAC));


thetaz = @(x) ((-ROy/2).*(x>=0).*(x.^2) + (FAt/2).*(x>=LOA).*((x-LOA).^2) + (RCy/2).*(x>=(LOA+LAC)).*((x-(LOA+LAC)).^2) ...
    -((FBr/2).*(x>=Ltot).*(x-Ltot).^2) + C3)/(E*I(d));
Yz = @(x) ((-ROy/6).*(x>=0).*(x.^3) + (FAt/6).*(x>=LOA).*((x-LOA).^3) + (RCy/6).*(x>=(LOA+LAC)).*((x-(LOA+LAC)).^3) ...
    -((FBr/6).*(x>=Ltot).*(x-Ltot).^3) + C3.*x + C4)/(E*I(d));

figure(5);
subplot(2,1,1);
plot(x,thetaz(x).*180/pi);
title('Slope along shaft in xy-plane [degrees]');
xlabel('Distance along the x axis (in)');
ylabel('Slope [degrees]');
subplot(2,1,2);
plot(x,Yz(x));
title('Deflection of shaft in xy-plane [in]');
xlabel('Distance along the x axis (in)')
ylabel('Deflection [in];')

step = 0:.001:Ltot;
XYslope = TrapInt(step,Mz(step)/(E*I(d)));
XZslope = TrapInt(step,My(step)/(E*I(d)));
XYdelta = TrapInt(step,XYslope);
XZdelta = TrapInt(step,XZslope);

%Subplot for slope and deflection for xy plane
figure(6);
subplot(2,1,1);
plot(step,XYslope.*(180/pi));
title('Slope along shaft in xy-plane [degrees]');
xlabel('Distance along the x axis (in)');
ylabel('Slope [degrees]');
subplot(2,1,2);
plot(step,XYdelta);
title('Deflection of shaft in xy-plane [in]');
xlabel('Distance along the x axis (in)')
ylabel('Deflection [in];')

figure(7);
subplot(2,1,1);
plot(step,XZslope.*(180/pi));
title('Slope along shaft in xz-plane [degrees]');
xlabel('Distance along the x axis (in)');
ylabel('Slope [degrees]');
subplot(2,1,2);
plot(step,XZdelta);
title('Deflection of shaft in xz-plane [in]');
xlabel('Distance along the x axis (in)')
ylabel('Deflection [in];')

newstep = 0:.001:30;
Iyop = find(newstep == 30);
c3XY = XYdelta(Iyop)/30;
c3XZ = XZdelta(Iyop)/30;

c2XY = atan(XYslope(Iyop)./(LOA+LAC));

c2XZ = atan(XZslope(Iyop)./(LOA+LAC)); 

correctedXYslope = XYslope - c2XY;
correctedXYdelta = XYdelta - c3XY.*step;
correctedXZslope = XZslope - c2XZ;
correctedXZdelta = XZdelta - c3XZ.*step;

magSlope = sqrt((XYslope).^2 + (XZslope).^2);
magDelta = sqrt((XYdelta).^2 + (XZdelta).^2);
cmagSlope = sqrt((XYslope - c2XY).^2 + (XZslope - c2XZ).^2);
cmagDelta = sqrt((XYdelta - c3XY.*step).^2 + (XZdelta - c3XZ.*step).^2);

figure(8);
subplot(2,1,1);
plot(step,correctedXYslope.*(180/pi)); hold on;
plot(step,XYslope.*(180/pi)); hold off;
title('Corrected Slope along shaft in xy-plane [degrees]');
xlabel('Distance along the x axis (in)');
ylabel('Slope [degrees]');
legend('Corrected','As Integrated');
subplot(2,1,2);
plot(step,correctedXYdelta); hold on;
plot(step,XYdelta); hold off;
title('Corrected Deflection of shaft in xy-plane [in]');
xlabel('Distance along the x axis (in)')
ylabel('Deflection [in];')
legend('Corrected','As Integrated');

figure(9);
subplot(2,1,1);
plot(step,correctedXZslope.*(180/pi)); hold on;
plot(step,XZslope.*(180/pi)); hold off;
title('Slope along shaft in xz-plane [degrees]');
xlabel('Distance along the x axis (in)');
ylabel('Slope [degrees]');
subplot(2,1,2);
plot(step,correctedXZdelta); hold on;
plot(step,XZdelta); hold off;
title('Corrected Deflection of shaft in xz-plane [in]');
xlabel('Distance along the x axis (in)')
ylabel('Deflection [in];');
legend('Corrected','As Integrated');

figure(10);
subplot(2,1,1);
plot(step,cmagSlope.*(180/pi)); hold on;
plot(step,magSlope.*(180/pi)); hold off;
title('Combined Magnitude of Slope [degrees]');
xlabel('Distance along the x axis (in)');
ylabel('Slope [degrees]');
subplot(2,1,2);
plot(step,cmagDelta); hold on;
plot(step,magDelta); hold off;
title('Combined Magnitude of Deflection [in]');
xlabel('Distance along the x axis (in)')
ylabel('Deflection [in];');
legend('Corrected','As Integrated');

%% Von Mises for cross section at bearing C
principals = ['principalA';'principalB';'principalC';'principalD'];
Von = [];
%eval(principals(1,:));
for i = 1:4
    
    next = getVonMises(eval(principals(i,:)));
    Von = [Von,next];
end
%
von1 = getVonMises(principalA);
von2 = getVonMises(principalB);
von3 = getVonMises(principalC);
von4 = getVonMises(principalD);

%% Step 8
% Determine the new minimum shaft diameter if the maximum slope at the two 
% bearing locations has to be less than 0.07 degrees.  (This ensures that the bearings do not bind)
rad_to_deg = 180/pi;
% general formula for moment of inertia for circular sross section
J = @(d) (pi/32)*d.^4; % area polar moment of inertia in inches ^4
I = @(d) (pi/64)*d.^4; % second area moment of inertia in inches ^4
d = 1.2; % initial diameter
locationOfInterest = find(step==30);
check =  cmagSlope(locationOfInterest).*(180/pi);
while check <= 0.07
    step = 0:.001:Ltot;
    XYslope = TrapInt(step,Mz(step)/(E*I(d)));
    XZslope = TrapInt(step,My(step)/(E*I(d)));
    XYdelta = TrapInt(step,XYslope);
    XZdelta = TrapInt(step,XZslope);
    newstep = 0:.001:30;
    Iyop = find(newstep == 30);
    c3XY = XYdelta(Iyop)/30;
    c3XZ = XZdelta(Iyop)/30;

    c2XY = atan(XYslope(Iyop)./(LOA+LAC));

    c2XZ = atan(XZslope(Iyop)./(LOA+LAC)); 

    correctedXYslope = XYslope - c2XY;
    correctedXYdelta = XYdelta - c3XY.*step;
    correctedXZslope = XZslope - c2XZ;
    correctedXZdelta = XZdelta - c3XZ.*step;

    magSlope = sqrt((XYslope).^2 + (XZslope).^2);
    magDelta = sqrt((XYdelta).^2 + (XZdelta).^2);
    cmagSlope = sqrt((XYslope - c2XY).^2 + (XZslope - c2XZ).^2);
    cmagDelta = sqrt((XYdelta - c3XY.*step).^2 + (XZdelta - c3XZ.*step).^2);

    
    
    d = d - 0.01;
    check = cmagSlope(locationOfInterest).*(180/pi)
end
%%
% some decisions have to be made on the design above
dOA = 1.2; % design gave
dAC = 1.2; % design gave 
dCB = 1.2; % design gave 

% find the polar moment of inertia all in inche ^4
JOA = J(dOA);
JAC = J(dAC);
JCB = J(dCB);
Ta = 0;
Tm = T;
% Angular deflection is just added so let us do that
Tmax = Ta+Tm; % the maximum torque seen

theta = (Tmax/G)*(Lac/JAC+Lcb/JCB); % deflection in radians
theta_deg = theta*rad_to_deg;

fprintf(1,'Deflection due to Torsion = %g deg\n',theta_deg);

%%
function sigma2 = getsigma2(principal)
    for i=1:3
        if max(principal) ~= principal(i) && min(principal) ~= principal(i)
            sigma2 = principal(i);
        end
    end
end

function VonMises = getVonMises(principal)
    %fprintf('Finding VonMises for: %\n',principal);
    sig1 = min(principal);
    sig2 = getsigma2(principal);
    sig3 = max(principal);
    VonMises = sqrt(.5*((sig1-sig2)^2+(sig2-sig3)^2+(sig3-sig1)^2));
end

function area = TrapInt(x,y)
% It takes the x and y values of a function and finds the area under the
% curve using the Trapeziodal rule
% input the steps of the independent variable x
% the values of the function evaluated at x
% Remember formula for area of trapezium 
% area = 0.5*sum(parallel sides)/height
% here the parallel sides are the value of the function
% and height is the step size !!

area(1) = 0;    % initialize summing variable
for i = 1:1:(length(x)-1)
    area(i+1) = area(i) + (1/2)*(y(i) + y(i+1))*(x(i+1) - x(i));
end
end