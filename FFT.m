clc; clear all; 
path = '/home/scott/GIT/PSU_ME_Stuff/ME345W/'; % Set the path to data files

%path = 'Enter Your Path to Data'
cd(path); % change the working directory
dir_list1 = dir('*.*'); % create a list of files in data folder 
file = 'FFT_homework.xlsx';
data = importdata(file);
freq = (0:1:63);
t = data.data(:,1);
Yt = data.data(:,2);
Yfft = fft(Yt);
k = (0:length(Yt)-1);
figure(1);
plot(freq,abs(Yfft));
figure(2);
plot(t,Yt);