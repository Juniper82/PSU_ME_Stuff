
function ICE4()
%% ICE 4 (the first part till 303 is the same as ICE 2)
% Shaft Design for Steady Torsion and Fully Reversed Bending
% Revision: Spring 2018
% This computes the defection of the shaft using singularity functions
% 
clear; clc; close all;
%sets the default line width in plots.  Needs to be done
%only once, but it is done here so we are sure.
set(0,'DefaultLineLineWidth',2);

%% some conversions defined
hp_to_lb_in_per_sec = 6600; % 1 hp = 6600 lb*in/sec
rpm_to_radper_sec = 2*pi/60; % 1 rmp = 2*pi/60 rads/sec
deg_to_rad = (2*pi / 360); % conversion from deg to radians
rad_to_deg = 360/(2*pi); % conversion from rad to deg
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
% Gear Pressure Angle
PHIdeg = 20;    % deg
PHI = PHIdeg * deg_to_rad   % rad

%% Solution
%
%% 1 Find the toruqe based on input power
T = P*hp_to_lb_in_per_sec/(OMEGA * rpm_to_radper_sec)     % lb in
% In other problems these may be different and will be needed to computed
% individually for the alternate and mean part
Ta = T; %In this problem as the alternating component
Tm = T; % There is also a mean torque

%% 2 Find the belt force and total belt force.  Again take care of the alternating and mean part
%Using the above Torque and belt tension ratio as 5
RHO = 5     % Tension ratio - no units
Fna = Ta/rs   % lb
F1a = Fna / (1 - (1/RHO)) % lb
F2a = F1a - Fna    % lb
Fsa = F1a + F2a    % lb

%similrly for the mean part
Fnm = Tm/rs   % lb
F1m = Fnm / (1 - (1/RHO)) % lb
F2m = F1m - Fnm    % lb
Fsm = F1m + F2m    % lb

%% 3 Find the force due to the gear
% Again find the alternating and mean part
Fgta = Ta / rg    % lb
Fgra = Fgta * tan(PHI)    % lb

Fgtm = Tm / rg    % lb
Fgrm = Fgtm * tan(PHI)    % lb
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
z = 0:0.01:8;    % in
% alternating component
qxza = @(z) R1xa .* (z==0) + Fgra .* (z==p) + R2xa .*(z==b) + Fsa .* (z==q);
Vxza = @(z) R1xa .* (z>=0).* z.^0 + Fgra .* (z>=p).*(z-p).^0 + ...
    R2xa .*(z>=b).* (z-b).^0 + Fsa .* (z>=q).*(z-q).^0;
Mya = @(z) R1xa .* (z>=0) .* z + Fgra .* (z>=p) .*...
    (z-p) + R2xa .* (z>=b) .* (z-b) + Fsa .* (z>=q) .* (z-q);
% mean component
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

subplot(3,1,2);
plot(z,Vxza(z),z,Vxzm(z),z,Vxza(z)+Vxzm(a));
grid on;
xlabel('Distance along the z axis (in)')
ylabel('Shear Force in x-z plane');
title('Shear Force in x-z plane');

subplot(3,1,3);
plot(z,Mya(z),z,Mym(z),z,Mya(z)+Mym(z));
xlabel('Distance along the z axis (in)')
ylabel('Moment in x-z plane (lb in)')
grid on
title('Moment in the x-z plane');

figure(myfig+1);
myfig=myfig+1;
Tm_arr=@(z) 0*(z>=0)+Tm.*(z>a)-Tm.*(z>c);
Ta_arr=@(z) 0*(z>=0)+Ta.*(z>a)-Ta.*(z>c);
plot(z,Tm_arr(z),z,Ta_arr(z),z,Ta_arr(z)+Tm_arr(z));
xlabel('Distance along the z axis (in)')
ylabel('Torsion in shaft(lb in)')
grid on
title('Torsion in shafts plane');

% similarily compute the forces in the yz place
% alternating component
qyza = @(z) R1ya .* (z==0) - Fgta .* (z==p) + R2ya .* (z==b);
Vyza = @(z) R1ya .*(z>=0).* z.^0 - Fgta .*(z>=p).* (z-p).^0 + ...
    R2ya .*(z>=b).* (z-b).^0;
Mxa = @(z) R1ya .* (z>=0) .* z - Fgta .* (z>=p) .*(z-p) + ...
    R2ya .* (z>=b) .* (z-b);

% mean component
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

subplot(3,1,2);
plot(z,Vyza(z),z,Vyzm(z),z,Vyza(z)+Vyzm(z));
grid on;
xlabel('Distance along the z axis (in)')
ylabel('Shear Force in y-z plane');
title('Shear Force in y-z plane');

subplot(3,1,3);
plot(z,Mxa(z),z,Mxm(z),z,Mxa(z)+Mxm(z))
xlabel('Distance along the z axis (in)')
ylabel('Moment in y-z plane (lb in)')
grid on
title('Moment in the y-z plane');

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

subplot(3,1,2);
Va = @(z)(Vxza(z).^2 + Vyza(z).^2).^(1/2);
Vm = @(z)(Vxzm(z).^2 + Vyzm(z).^2).^(1/2);
plot(z,Va(z),z,Vm(z),z,Va(z)+Vm(z));
grid on;
xlabel('Distance along the z axis (in)')
ylabel('Shear Force');
title('Magnitude of Shear Force');

subplot(3,1,3);
Ma = @(z)(Mxa(z).^2 + Mya(z).^2).^(1/2);
Mm = @(z)(Mxm(z).^2 + Mym(z).^2).^(1/2);
plot(z,Ma(z),z,Mm(z),z,Ma(z)+Mm(z))
xlabel('Distance along the z axis (in)')
ylabel('Total Moment magnitude (lb in)')
grid on
title('Combined Moment magnitude');

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
% put z back to the full length of the shaft
z = 0:0.001:8;    % in
%% 7 Using a trial material obtain the corrections for the factigue loading

Sut = 68*1000;   % psi
Sy = 38*1000;    % psi
Sprme = 0.5 * Sut;   % psi
fprintf(1,'Uncorrected Endurance Strength = %g psi\n',Sprme);

Cload = 1;
Csize = @(d) 0.869 * d^-0.097;  % no units (assumes d [=] inches)
Csurf = 2.7 * (Sut/1000)^-0.265;        % no units (assumes Sut [=] ksi)
Ctemp = 1;
Creliab = 1;

Se = @(d) Cload * Csize(d) * Csurf *...
    Ctemp * Creliab * Sprme;

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
fprintf(1,'Problem we have stated really begins here.  Press Any Key...\n');
pause;
%% ICE 3 really begins here
%%%%%%%%%%%%%ICE 3 really starts here, but all the previous calculations
%%%%%%%%%%%%%%%are needed
% Material properties obtained from 
% http://matweb.com/search/DataSheet.aspx?MatGUID=53ad0dbdb3c143d89c38c47904e3f2ff&ckck=1
% and rounded up or down
%% Basic material properties
G = 11.7e6; % modulus of rigidity in psi
Ym = 29.9e6; % Young Modulus in psi
lAB = 1.5; %length in inches
lBC = 3.5; %length in inches
lCD = 1.5; %length in inches
lDE = 1.5; %length in inches
%% general formula for moment of inertia for circular sross section
J = @(d) (pi/32)*d.^4; % area polar moment of inertia in inches ^4
In = @(d) (pi/64)*d.^4; % second area moment of inertia in inches ^4

% some decisions have to be made on the design above
dAB = 0.75; % Nominal diaof material we choose
dBC = 0.6281; % design gave 
dCD = 0.6112; % design gave 
dDE = 0.5077; % should have 

% find the polar moment of inertia all in inche ^4
JAB = J(dAB);
JBC = J(dBC);
JCD = J(dCD);
JDE = J(dDE);
% Angular deflection is just added so let us do that
Tmax = Ta+Tm; % the maximum torque seen

theta = (Tmax/G)*(lAB/JAB+lBC/JBC+lCD/JCD); % deflection in radians
theta_deg = theta*rad_to_deg;

fprintf(1,'Deflection due to Torsion = %g deg\n',theta_deg);
fprintf(1,'Press Any Key...\n');
pause;
%%  ICE 4 begins here to find the slope and deflection 
% Slope and Delfection
% Now find the slope and deflection of the shaft - remember we have to 
% deal with this in two differenr places!!

%Define the length of the shaft
STEP = 0.001;
z1=[0:STEP:1.5];
z2=[1.5:STEP:5];
z3=[5:STEP:6.5];
z4=[6.5:STEP:8.0];

% find the moment of inertia for each section
InAB = In(dAB)
InBC = In(dBC)
InCD = In(dCD)
InDE = In(dDE)

% % Consider the total moment for this part
% M = @(z) Ma(z)+Mm(z);
% Consider the total moment for this part
Mx = @(z) Mxa(z)+Mxm(z); % combine the alternate and mean moment in x
My = @(z) Mya(z)+Mym(z); % combine the alternate and mean moment in y
% Ma = @(z) sqr(Mxa(z).^2+Mya(z).^2);
% Mm = @(z) sqrt(Mxm(z).^2+Mym(z).^2);
% MM = @(z) Ma(z)+Mm(z);
% Now plot the M/EI function
figure(myfig);
myfig=myfig+1;
plot([z1,z2,z3,z4],[Mx(z1)/(Ym*InAB),Mx(z2)/(Ym*InBC),...
    Mx(z3)/(Ym*InCD),Mx(z4)/(Ym*InDE)]);
grid on;
disp('M/EI function');
xlabel('Distance (in)');
ylabel('M/EI function');
% Now find the slope in the yz plane
s1x = TrapInt(z1,(Mx(z1)/(Ym*InAB)));
s2x = s1x(end)+TrapInt(z2,(Mx(z2)/(Ym*InBC)));
s3x = s2x(end)+TrapInt(z3,(Mx(z3)/(Ym*InCD)));
s4x = s3x(end)+TrapInt(z4,(Mx(z4)/(Ym*InDE)));
% Now find the slope in the xz plane
s1y = TrapInt(z1,(My(z1)/(Ym*InAB)));
s2y = s1y(end)+TrapInt(z2,(My(z2)/(Ym*InBC)));
s3y = s2y(end)+TrapInt(z3,(My(z3)/(Ym*InCD)));
s4y = s3y(end)+TrapInt(z4,(My(z4)/(Ym*InDE)));

% s1 = TrapInt(z1,(MM(z1)/(Ym*InAB)));
% s2 = s1(end)+TrapInt(z2,(MM(z2)/(Ym*InBC)));
% s3 = s2(end)+TrapInt(z3,(MM(z3)/(Ym*InCD)));
% s4 = s3(end)+TrapInt(z4,(MM(z4)/(Ym*InDE)));

% figure(100);
% ss = sqrt(([s1x,s2x,s3x,s4x]).^2+([s1y,s2y,s3y,s4y]).^2);
% plot([z1,z2,z3,z4],ss,[z1,z2,z3,z4],[s1,s2,s3,s4]);

% Now find the deflection in the yz plane
d1x = TrapInt(z1,s1x);
d2x = d1x(end)+TrapInt(z2,s2x);
d3x = d2x(end)+TrapInt(z3,s3x);
d4x = d3x(end)+TrapInt(z4,s4x);
% Now find the deflection in the xz plane
d1y = TrapInt(z1,s1y);
d2y = d1y(end)+TrapInt(z2,s2y);
d3y = d2y(end)+TrapInt(z3,s3y);
d4y = d3y(end)+TrapInt(z4,s4y);

% Now we need to be sure that deflection at z=5 should be zero, but it will
% % not and so we apply a correction
I = find(z2==5);

c3sx = d2x(I)/5;
c3sy = d2y(I)/5;

c2sx = atan(s2x(I))/5;
c2sy = atan(s2y(I))/5;

%total magnitude of slope and magnitude of deflection after all correction
ss = sqrt(([s1x,s2x,s3x,s4x]-c2sx).^2+([s1y,s2y,s3y,s4y]-c2sy).^2);
ddd = sqrt(([d1x,d2x,d3x,d4x]-c3sx*[z1,z2,z3,z4]).^2+...
    ([d1y,d2y,d3y,d4y]-c3sy*[z1,z2,z3,z4]).^2);

% now plot everything

figure(myfig);
myfig = myfig+1;
plot([z1,z2,z3,z4],[s1x,s2x,s3x,s4x]);hold on;
plot([z1,z2,z3,z4],[s1x,s2x,s3x,s4x]-c2sx);hold off;
xlabel('Distance (in)');
ylabel('Slope (radians)');
grid on;
legend('As Integrated','Corrected');
title('y slope');
%%
fprintf(1,'Press Any Key...\n');
pause;
%%
figure(myfig);
myfig = myfig+1;
plot([z1,z2,z3,z4],[s1y,s2y,s3y,s4y]);hold on
plot([z1,z2,z3,z4],[s1y,s2y,s3y,s4y]-c2sy);hold off;
xlabel('Distance (in)');
ylabel('Slope (radians)');
grid on;
title('x slope');
legend('As Integrated','Corrected');
%%
fprintf(1,'Press Any Key...\n');
pause;
%%

figure(myfig);
myfig = myfig+1;
plot([z1,z2,z3,z4],ss);
xlabel('Distance (in)');
ylabel('Slope (radians)');
title('Magnitude of slope');
grid on;
%%
fprintf(1,'Press Any Key...\n');
pause;
%%

figure(myfig);
myfig = myfig+1;
plot([z1,z2,z3,z4],[d1x,d2x,d3x,d4x]);hold on
plot([z1,z2,z3,z4],[d1x,d2x,d3x,d4x]-c3sx*[z1,z2,z3,z4]);
xlabel('Distance (in)');
ylabel('deflection (in)');
grid on;
title('y deflection');
legend('As Integrated','Corrected');
%%
fprintf(1,'Press Any Key...\n');
pause;
%%

figure(myfig);
myfig = myfig+1;
plot([z1,z2,z3,z4],[d1y,d2y,d3y,d4y]);hold on;
plot([z1,z2,z3,z4],[d1y-c3sy*z1,d2y-c3sy*z2,d3y-c3sy*z3,d4y-c3sy*z4]);
xlabel('Distance (in)');
ylabel('deflection (in)');
legend('As Integrated','Corrected');
grid on;
title('x deflection');
hold off;
%%
fprintf(1,'Press Any Key...\n');
pause;
%%

% plot the total deflection
figure(myfig);
myfig = myfig+1;
plot([z1,z2,z3,z4],ddd);
grid on;
xlabel('Distance (in)');
ylabel('deflection (in)');
title('Magnitude of Deflection');
legend('As Integrated','Corrected');
% put up texts
text(0,0.0001,'A');
text(p,0.0001,'B');
text(b,0.0001,'C');
text(c,0.0001,'D');
%%
fprintf(1,'Press Any Key...\n');
pause;
%%

% Now we need to be sure that deflection at z=5 should be zero, but it will
% % not and so we apply a correction
I = find(z2==5);
c3sx = d2x(I)/5;
c3sy = d2y(I)/5;

c2sx = s2x(I)/5;
c2sy = s2y(I)/5;
ss(find(z>b,1))
fprintf(1, 'slope @C = %g rad, slope @B = %g rad, slope @D = %g rad\n',...
    ss(find(b<z,1)),ss(find(p<z,1)),ss(find(c<z,1)));
fprintf(1, 'delta @C = %g in, delta @B = %g in, delta @D = %g in\n',...
    ddd(find(b<z,1)),ddd(find(p<z,1)),ddd(find(c<z,1)));
fprintf(1,'================================================================\n');
%%
fprintf(1,'Press Any Key...\n');
pause;
%%

%%
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