function [spikeMat, tVec] = Coding(tSim_us, InMatrix_mA)
% Poisson coding 
load('param.mat');
tVec = 0:dt:tSim_us-dt;
nBins = floor(tSim_us/dt);
nTrials=size(InMatrix_mA,2);
InMatrix_mA=InMatrix_mA*100;
for i=1:nTrials
    spikeMat(i,:) = rand(1, nBins) < InMatrix_mA(i)*0.1*dt;
end