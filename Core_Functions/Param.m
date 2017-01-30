clc
clear

%% Coding parameters
dt = 1/10; % Spike accuracy

%% SRM parameters
vthr=1;% SRM : neuron threshold 
tav=0.05;% SRM : tav
amp=1;%SRM : alpha function amplitude 
nu=0;%SRM : Refractory priode amplitude

%% STDP parameters
wmax=1;
eta=0.08;
A_LTP = eta*0.4;
A_LTD = eta*0.1;
tav_LTP = 0.04;
tav_LTD = 0.01;

save('F:\PHD sience\Simulation\COM\Functions\param.mat');