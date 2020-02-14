clc; clear all;
p = 2.0;
b = 5.0;
q = 6.75;

R1x = 7.65;
R2x = -52.87;
R1y = 14.64;
R2y = 9.74;
Fgr = 8.78;
Fgt = 24.36;
Fs = 36.54;

%%
z = 0:0.01:6.75;    % "

% Singularity functions for x-z plane
qxz = @(z) R1x .* (z==0) + Fgr .* (z==p) + R2x.*(z==b) + Fs .* (z==q);
Vxz = @(z) R1x .* (z>=0).*z.^0 + Fgr.*(z>=p).*(z-p).^0 + R2x .* (z>=b).*(z-b).^0 ...
    + Fs.*(z>=q).*(z-q);
Mxz = @(z) R1x .* (z>=0).*z + Fgr.*(z>=p).*(z-p) + R2x .* (z>=b).*(z-b) ...
    + Fs.*(z>=q).*(z-q);

% Singularity functions for y-z plane
qyz = @(z) R1y .* ( z==0) - Fgt .*(z==p) + R2y .* (z==b);
Vyz = @(z) R1y .*(z>=0).*z.^0 - Fgt .*(z>=p).*(z-p).^0 + R2y .* (z>=b).*(z-b).^0;
Myz = @(z) R1y .*(z>=0).*z - Fgt .*(z>=p).*(z-p) + R2y .* (z>=b).*(z-b);

%%
figure(1);

subplot(3,1,1);
plot(z,qxz(z)); grid on;
xlabel('Distance along the z axis (in)')
ylabel('Force in x-z plane');
title('Force in x-z plane');
%
subplot(3,1,2);
plot(z,Vxz(z));
grid on;
xlabel('Distance along the z axis (in)')
ylabel('Shear Force in x-z plane (lb)');
title('Shear Force in x-z plane');
%
subplot(3,1,3);
plot(z,Mxz(z));
xlabel('Distance along the z axis (in)')
ylabel('Moment in x-z plane (lb.in)')
grid on
title('Moment in the x-z plane');

figure(2);

subplot(3,1,1);
plot(z,qyz(z)); grid on;
xlabel('Distance along the z axis (in)')
ylabel('Force in y-z plane');
title('Force in y-z plane');
%
subplot(3,1,2);
plot(z,Vyz(z));
grid on;
xlabel('Distance along the z axis (in)')
ylabel('Shear Force in y-z plane (lb)');
title('Shear Force in y-z plane');
%
subplot(3,1,3);
plot(z,Myz(z));
xlabel('Distance along the z axis (in)')
ylabel('Moment in y-z plane (lb.in)')
grid on
title('Moment in the y-z plane');
A = 0;
B = 2;
C = 5;
D = 6.5;

M_A = sqrt(Mxz(A)^2 + Myz(A)^2)
M_B = sqrt(Mxz(B)^2 + Myz(B)^2)
M_C = sqrt(Mxz(C)^2 + Myz(C)^2)
M_D = sqrt(Mxz(D)^2 + Myz(D)^2)