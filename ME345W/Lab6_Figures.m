% Lab 6 ME 345W 
clc; clear all;

% Change working directory to directory containing lab 2 data and read it
% in as a structure

path = '/home/scott/Documents/ME345W/Labs/Lab6'; % Set the path to data files

%path = 'Enter Your Path to Data'
cd(path); % change the working directory
dir_list1 = dir('*.*'); % create a list of files in data folder 
file = 'Lab6.csv';
data = importdata(file);

%%
Temp = data.data(:,1) + 273.15;
Resistances = data.data(:,2).*10^3;

A = [ones(length(Resistances),1), log(Resistances), log(Resistances).^3];
b = 1./Temp;
% Below are two simple approaches to solving a linear system in Matlab
x = linsolve(A,b);
x1 = A\b;
% Matlab is giving me a warning for Dr. Litwhiler's code
X = inv(A'*A)*A'*b;
% All 3 produce methods to solve the linear system produce the same results
C1 = X(1,1);
C2 = X(2,1);
C3 = X(3,1);

fprintf('C1 = %.4e\nC2 = %.4e\nC3 = %.4e\n',C1,C2,C3)

figure(1);
plot(Temp, Resistances,'*k','MarkerSize',8)
grid on;
ylabel('\bf Resistance [\Omega]')
xlabel('\bf Temperature [K]')
title({'Analysis of an Omega 44006 Thermistor','Resistance vs Temperature'})
export_fig /home/scott/Documents/ME345W/Labs/Lab6/Figure1.png -transparent
%%
file = 'Lab6_set2.csv';
dat = importdata(file);
R = dat.data(:,1)*1000;

T1 = dat.data(:,2);
T2 = dat.data(:,3);
T3 = dat.data(:,4);
figure(2);
plot(R,T1,'.r',R,T2,'^k',R,T3,'*k','MarkerSize',8);
title({'Comparing Expected and Experimental Values of 1/T',''});
ylabel('\bf 1/T [1/K]');
xlabel('\bf Resistance [\Omega]');
legend('Measured 1/T','Data Sheet Coefficients Eqn(1)','Experimental Coefficients Eqn(1)','Location','NorthWest');
grid on;
export_fig /home/scott/Documents/ME345W/Labs/Lab6/Figure2.png -transparent