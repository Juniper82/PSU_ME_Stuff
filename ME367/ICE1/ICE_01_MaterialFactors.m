%% ICE 1
% Shaft Design for Steady Torsion and Fully Reversed Bending
% Revision: Spring 2019
% Now select the material and find all the factors etc.
%% 7 Using a trial material obtain the corrections for the factigue loading

% For SAE 1030 hot rolled from page 1024, table A-9
Sut = 68*1000;   % psi
Sy = 38*1000;    % psi
Sprme = 0.5 * Sut;   % psi

Cload = 1;% eq 6.7a, pp362
Csize = @(d) 0.869 * d^-0.097;  % eq 6.7 b p 363, no units (assumes d [=] inches)
% Important this equation uses the Sut in ksi
Csurf = 2.7 * (Sut/1000)^-0.265; % eq 6.7e p 365 OR table 6.3  no units (assumes Sut [=] ksi)
Ctemp = 1; % eq 6.7f pp 367
Creliab = 1; % table 6.4, p367

% corrected endurance limit
Se = @(d) Cload * Csize(d) * Csurf * Ctemp * Creliab * Sprme;

%
%% 8 : Calculate the notch sensitivity
% Notch sensitivity for bending 
a = 0.100^2     % in
r = 0.010;      % in
q = 1 / (1 + sqrt(a/r)) % eq 6.35, pp 376

% notch sensitivity for torsion

a = 0.075^2     % in
r = 0.010;      % in
qs = 1 / (1 + sqrt(a/r)) % eq 6.35, pp 376

%% 9 Obtain the fatigue stress concentration factor for each section for tension and shear and keyway
% These are given values in this problem.  In general you do know diameter
% and radius to determine this.  So you start with some guess values
Ktb = 3.5;
Kts = 2.0;
Ktk = 4.0;
disp('Fatigue stress concentration for section C :');
Kfb = 1 + q * (Ktb - 1);
Kfs = 1 + qs * (Kts - 1);

Kfsm = Kfs;
fprintf(1,'Kfb = %g, Kfs = %g, Kfsm = %g\n',Kfb, Kfs, Kfsm);
