%% ICE 2
% Shaft Design for Steady Torsion and Fully Reversed Bending
% Revision: Spring 2019
%clc;clear all;
%% 7 Using a trial material obtain the corrections for the factigue loading
% For SAE 1030 hot rolled from page 1024, table A-9
Sut = 68*1000;   % psi
Sy = 38*1000;    % psi
Sprme = 0.5 * Sut;   % psi
fprintf(1,'Uncorrected Endurance Strength = %g psi\n',Sprme);

Cload = 1;
Csize = @(d) 0.869 * d^-0.097;  % no units (assumes d [=] inches)
%SurfaceFinish = "CR"; % cold rolled OR Machines
SurfaceFinish = "HR"; % hot rolled
if strcmp(SurfaceFinish,"CR")
    Asurf = 2.7; %No units
    bsurf = -0.265; % No units
end
if strcmp(SurfaceFinish,"HR")
    Asurf = 14.4; %No units
    bsurf = -0.718; % No units
end;

Csurf = Asurf * (Sut/1000)^bsurf; % eq 6.7e p 365 OR table 6.3  no units (assumes Sut [=] ksi)
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

