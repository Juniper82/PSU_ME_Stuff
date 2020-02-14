%% ICE 2
% Shaft Design for Steady Torsion and Fully Reversed Bending
% Revision: Spring 2020
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
rg = 0.5*dg     % in
ds = 6.00;      % in
rs = 0.5*ds     % in
fprintf(1,'Gear Radius = %g in, Sheave Radisu = %g in \n',rg,rs);
% Gear Pressure Angle
PHIdeg = 20;    % deg
PHI = PHIdeg * deg_to_rad;   % rad
fprintf(1,'PHI = %g radians\n',PHI);
%% Solution
%% 1 Find the toruqe based on input power
T = P*hp_to_lb_in_per_sec/(OMEGA * rpm_to_radper_sec);     % lb in
fprintf(1,'T = %g lb*in\n',T);
%% Solution
%
%% 1 Find the toruqe based on input power
T = P*hp_to_lb_in_per_sec/(OMEGA * rpm_to_radper_sec)     % lb in
% In other problems these may be different and will be needed to be computed
% individually for the alternate and mean part.  These could also be
% functions and you will have to suitably handle them
Ta = T; %In this problem as the alternating component
Tm = T; % There is also a mean torque

%% 2 Find the belt force and total belt force.  Again take care of the alternating and mean part
%Using the above Torque and belt tension ratio as 5
RHO = 5;     % Tension ratio - no units
Fna = Ta/rs;   % lb
F1a = Fna / (1 - (1/RHO)); % lb
F2a = F1a - Fna;    % lb
Fsa = F1a + F2a;    % lb

fprintf(1,'Fna = %g lb, F1a = %g lb\n',Fna,F1a);
fprintf(1,'F2a = %g lb, Fsa = %g lb\n',F2a,Fsa);

%similrly for the mean part
Fnm = Tm/rs;   % lb
F1m = Fnm / (1 - (1/RHO)); % lb
F2m = F1m - Fnm;    % lb
Fsm = F1m + F2m;    % lb
fprintf(1,'Fnm = %g lb, F1m = %g lb\n',Fnm,F1m);
fprintf(1,'F2m = %g lb, Fsm = %g lb\n',F2m,Fsm);

%% 3 Find the force due to the gear
% Again find the alternating and mean part
Fgta = Ta / rg    % lb
Fgra = Fgta * tan(PHI)    % lb

Fgtm = Tm / rg    % lb
Fgrm = Fgtm * tan(PHI)    % lb
fprintf(1,'Fgta = %g lb, Fgra = %g lb\n',Fgta,Fgra);
fprintf(1,'Fgtm = %g lb, Fgrm = %g lb\n',Fgtm,Fgrm);
%
%% 4 Assuming gear and sheave force to be concentrated at their center
% remember that these are distributed forces in reality
% Below reactions are based on the FBD drawn
% 
%The alternating part
R2xa = (1/b)*(-Fgra*p - Fsa*q)  % lb
R1xa = -Fgra - R2xa - Fsa   % lb
R2ya = (1/b)*(Fgta*p) % lb
R1ya = Fgta - R2ya     % lb

% the mean part
R2xm = (1/b)*(-Fgrm*p - Fsm*q)  % lb
R1xm = -Fgrm - R2xm - Fsm   % lb
R2ym = (1/b)*(Fgtm*p) % lb
R1ym = Fgtm - R2ym     % lb

%% 5 Discontunity functions and obtaining the shear force and bending moments
% using the discontinuity function (as shown in the handout)
% we will then plot the loading function, shear force and bending moment in
% each place and then finally the total magnitude

% define the load function
z = 0:0.01:7;    % in
qxza = @(z) R1xa .* (z==0) + Fgra .* (z==p) + R2xa .*(z==b) + Fsa .* (z==q);
Vxza = @(z) R1xa .* (z>=0).* z.^0 + Fgra .* (z>=p).*(z-p).^0 + ...
    R2xa .*(z>=b).* (z-b).^0 + Fsa .* (z>=q).*(z-q).^0;
Mya = @(z) R1xa .* (z>=0) .* z + Fgra .* (z>=p) .*...
    (z-p) + R2xa .* (z>=b) .* (z-b) + Fsa .* (z>=q) .* (z-q);

qxzm = @(z) R1xm .* (z==0) + Fgrm .* (z==p) + R2xm .*(z==b) + Fsm .* (z==q);
Vxzm = @(z) R1xm .* (z>=0).* z.^0 + Fgrm .* (z>=p).*(z-p).^0 + ...
    R2xm .*(z>=b).* (z-b).^0 + Fsm .* (z>=q).*(z-q).^0;
Mym = @(z) R1xm .* (z>=0) .* z + Fgrm .* (z>=p) .*...
    (z-p) + R2xm .* (z>=b) .* (z-b) + Fsm .* (z>=q) .* (z-q);


myfig=1;
figure(myfig);
myfig=myfig+1;
subplot(3,1,1);
plot(z,qxza(z),z,qxzm(z),z,qxzm(z)+qxza(z)); grid on;
xlabel('Distance along the z axis (in)')
ylabel('Force in x-z plane');
title('Force in x-z plane');
legend('Alt.','Mean','Total');

subplot(3,1,2);
plot(z,Vxza(z),z,Vxzm(z),z,Vxza(z)+Vxzm(a));
grid on;
xlabel('Distance along the z axis (in)')
ylabel('Shear Force in x-z plane');
title('Shear Force in x-z plane');
legend('Alt.','Mean','Total');

subplot(3,1,3);
plot(z,Mya(z),z,Mym(z),z,Mya(z)+Mym(z));
xlabel('Distance along the z axis (in)')
ylabel('Moment in x-z plane (lb in)')
grid on
title('Moment in the x-z plane');
legend('Alt.','Mean','Total');

figure(myfig);
myfig=myfig+1;
Tm_arr=@(z) 0*(z>=0)+Tm.*(z>a)-Tm.*(z>c);
Ta_arr=@(z) 0*(z>=0)+Ta.*(z>a)-Ta.*(z>c);
plot(z,Tm_arr(z),z,Ta_arr(z),z,Ta_arr(z)+Tm_arr(z));
xlabel('Distance along the z axis (in)')
ylabel('Torsion in shaft(lb in)')
grid on
title('Torsion in shafts plane');
legend('Alt.','Mean','Total');

% similarily compute the forces in the yz place
qyza = @(z) R1ya .* (z==0) - Fgta .* (z==p) + R2ya .* (z==b);
Vyza = @(z) R1ya .*(z>=0).* z.^0 - Fgta .*(z>=p).* (z-p).^0 + ...
    R2ya .*(z>=b).* (z-b).^0;
Mxa = @(z) R1ya .* (z>=0) .* z - Fgta .* (z>=p) .*(z-p) + ...
    R2ya .* (z>=b) .* (z-b);


qyzm = @(z) R1ym .* (z==0) - Fgtm .* (z==p) + R2ym .* (z==b);
Vyzm = @(z) R1ym .*(z>=0).* z.^0 - Fgtm .*(z>=p).* (z-p).^0 + ...
    R2ym .*(z>=b).* (z-b).^0;
Mxm = @(z) R1ym .* (z>=0) .* z - Fgtm .* (z>=p) .*(z-p) + ...
    R2ym .* (z>=b) .* (z-b);

figure(myfig)
myfig=myfig+1;
subplot(3,1,1);
plot(z,qyza(z),z,qyzm(z),z,qyzm(z)+qyza(z)); grid on;
xlabel('Distance along the z axis (in)')
ylabel('Force in y-z plane');
title('Force in y-z plane');
legend('Alt.','Mean','Total');

subplot(3,1,2);
plot(z,Vyza(z),z,Vyzm(z),z,Vyza(z)+Vyzm(z));
grid on;
xlabel('Distance along the z axis (in)')
ylabel('Shear Force in y-z plane');
title('Shear Force in y-z plane');
legend('Alt.','Mean','Total');

subplot(3,1,3);
plot(z,Mxa(z),z,Mxm(z),z,Mxa(z)+Mxm(z))
xlabel('Distance along the z axis (in)')
ylabel('Moment in y-z plane (lb in)')
grid on
title('Moment in the y-z plane');
legend('Alt.','Mean','Total');

figure(myfig);
myfig=myfig+1;
subplot(3,1,1);
Qa = @(z) ((qxza(z).^2 + qyza(z).^2).^(1/2));
Qm = @(z) ((qxzm(z).^2 + qyzm(z).^2).^(1/2));
plot(z,Qa(z),z,Qm(z),z,Qa(z)+Qm(z));
grid on;
xlabel('Distance along the z axis (in)')
ylabel('Force');
title('Magnitude of Force');
legend('Alt.','Mean','Total');

subplot(3,1,2);
Va = @(z)(Vxza(z).^2 + Vyza(z).^2).^(1/2);
Vm = @(z)(Vxzm(z).^2 + Vyzm(z).^2).^(1/2);
plot(z,Va(z),z,Vm(z),z,Va(z)+Vm(z));
grid on;
xlabel('Distance along the z axis (in)')
ylabel('Shear Force');
title('Magnitude of Shear Force');
legend('Alt.','Mean','Total');

subplot(3,1,3);
Ma = @(z)(Mxa(z).^2 + Mya(z).^2).^(1/2);
Mm = @(z)(Mxm(z).^2 + Mym(z).^2).^(1/2);
plot(z,Ma(z),z,Mm(z),z,Ma(z)+Mm(z))
xlabel('Distance along the z axis (in)')
ylabel('Total Moment magnitude (lb in)')
grid on
title('Combined Moment magnitude');
legend('Alt.','Mean','Total');

%%  6 obtain the bending moment at sections of interest
z = p         % in
MBa = Ma(z)   % lb in
MBm = Mm(z)   % lb in

z = b         % in
MCa = Ma(z)   % lb in
MCm = Mm(z)

z = c         % in
MDa = Ma(z)   % lb in
MDm = Mm(z)
%% 7 Using a trial material obtain the corrections for the factigue loading
% For SAE 1030 hot rolled from page 1024, table A-9
Sut = 68*1000;   % psi
Sy = 38*1000;    % psi
Sprme = 0.5 * Sut;   % psi
fprintf(1,'Uncorrected Endurance Strength = %g psi\n',Sprme);

Cload = 1;
Csize = @(d) 0.869 * d^-0.097;  % no units (assumes d [=] inches)
Csurf = 2.7 * (Sut/1000)^-0.265;        % no units (assumes Sut [=] ksi)
Ctemp = 1;
Creliab = 1;

Se = @(d) Cload * Csize(d) * Csurf * Ctemp * Creliab * Sprme;

%
%% 8 : Calculate the notch sensitivity
% Notch sensitivity for bending
a = 0.100^2;     % in
r = 0.010;      % in
q = 1 / (1 + sqrt(a/r));

% notch sensitivity for torsion

a = 0.075^2;     % in
r = 0.010;      % in
qs = 1 / (1 + sqrt(a/r));

%% 9 Obtain the fatigue stress concentration factor for each section for tension and shear
Ktb = 3.5;  %bending at a step
Kts = 2.0;  %torsion at a step
Ktk = 4.0;  %at a keyway

% Obtain the fatigue stress concentration factor based on the notch
% sensitivity and static stress concentration factors
Kfb = 1 + q * (Ktb - 1);
Kfs = 1 + qs * (Kts - 1);
fprintf(1,'Fatigue Stress Conc. @C Bending = %g, Torsion = %g\n',Kfb, Kfs);
Kfsm = Kfs;
Kfbm = Kfb;

%
%% 10 Obtain the shaft diameter based on trial guess at various locations
% safety factor
Nf = 2.5;
fprintf(1,'Factor of Safety = %f\n',Nf);

disp(' ');
disp('The required shaft diameter at point C:')
FUNC = @(d) d-(((32*Nf)/pi)*...
    (((sqrt((Kfb*MCa)^2 + (3/4)*(Kfs*Ta)^2))/Se(d))+...
    ((sqrt((Kfbm*MCm)^2 + (3/4)*(Kfsm*Tm)^2))/Sut)))^(1/3);
dc = fzero(FUNC,0.5);
disp(['Dia @ C = ', num2str(dc),' inches']);
% disp(['Csize(d) =', num2str(Csize(dc))]);    % no units
% disp(['Se(d) =', num2str(Se(dc)), 'psi']);    % psi

disp(' ');
disp('The required shaft diameter at point B:');
Kfb = 1 + q * (Ktk - 1);
Kfs = 1 + qs * (Ktk - 1);

Kfsm = Kfs;
Kfbm = Kfb;
fprintf(1,'Fatigue Stress Conc. @B (keyway) Bending = %g, Torsion = %g\n',Kfb, Kfs);
FUNC = @(d) d-(((32*Nf)/pi)*...
    (((sqrt((Kfb*MBa)^2 + (3/4)*(Kfs*Ta)^2))/Se(d))+...
    ((sqrt((Kfbm*MBm)^2 + (3/4)*(Kfsm*Tm)^2))/Sut)))^(1/3);
db = fzero(FUNC,0.5);
disp(['Dia @ B = ', num2str(db),' inches']); 
% disp(['Csize(d) =', num2str(Csize(db))])    % no units
% disp(['Se(d) =', num2str(Se(db)), 'psi'])    % psi

disp(' ');
disp('The required shaft diameter at point D:');

Kfb = 1 + q * (Ktb - 1);
Kfs = 1 + qs * (Kts - 1);

Kfsm = Kfs;
Kfbm = Kfb;
fprintf(1,'Fatigue Stress Conc. @D Bending = %g, Torsion = %g\n',Kfb, Kfs);
FUNC = @(d) d-(((32*Nf)/pi)*...
    (((sqrt((Kfb*MDa)^2 + (3/4)*(Kfs*Ta)^2))/Se(d))+...
    ((sqrt((Kfbm*MDm)^2 + (3/4)*(Kfsm*Tm)^2))/Sut)))^(1/3);
dd = fzero(FUNC,0.5);
disp(['Dia @ D = ', num2str(dd),' inches']);
% disp(['Csize(d) =', num2str(Csize(dd))])    % no units
% disp(['Se(d) =', num2str(Se(dd)), 'psi'])    % psi

fprintf(1,'================================================================\n');
fprintf(1,'Dia @ C = %g, Dia @ B = %g, Dia @ D = %g\n',dc,db,dd);
fprintf(1, 'Se @C = %g psi, Se @B = %g psi, Se @D = %g psi\n',Se(dc),Se(db),Se(dd));
fprintf(1,'================================================================\n');

