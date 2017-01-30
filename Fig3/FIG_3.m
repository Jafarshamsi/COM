% This code inclues :
%                   1. Generation of the dataset;
%                   2. Adding noise to the dataset;
%                   3. Coding based on Poisson distribution;
%                   4. STDP training of the WTA;
%                   5. Output of the WTA for Original data.
%                   6. Ploting the output rate of the neurons for each pattern
%                   7. Ploting the Recognition rate of the Original and noisy dataset.
% How to use:
%                   1. set the parameters
%                   2. RUN
%% CLEAR

clear
clc

%% parameters

n=6; % Number of patterns = Number of neurons
tT=20;% Training time
tR=500;% Retrival time
mu_weight_init=0.2;% Mu initial weight
sigma_weight_init=0.01;% Sigma initial weight

%% Generation of Dataset

 %Generation of Random Pattern (Linearly seperatable)
        l=14;
        sl=sqrt(l);
        for i=1:n
            data(i,1:l)=(0.8+0.1*rand(1,l)).*(-floor(-full(sprand(1,l,0.5))));% generate random dataset through Gaussian distribution
            data(i,l+1)=i;
        end
          
 
%% Coding 

for i=1:n
    Coding_Cnt=i
    [spikeMat_Coding1, tVec_Coding] = Coding(tT,data(i,1:l));
    spikeMat_Coding{i}(:,:)=spikeMat_Coding1;
end


%% Training

w_init_STDP=mu_weight_init + sigma_weight_init.*randn(1, l);
for i=1:n
    Training_Cnt=i
    [w_stdp1,tVec_stdp] = STDP(spikeMat_Coding{i},tVec_Coding,w_init_STDP);
    w_stdp{i}(:,:)=w_stdp1;
end

%% Processing Original data: Neurons outputs, Correct recognition (Crct), input rate (rate_in) and output rate (rate_out) of neurons
wstep=10;
Crct=zeros(1,wstep);
step=size(tVec_Coding,2)/10;
for s=1:wstep%select weights
    for i=1:n
        w_stdp_wta(i,:)=w_stdp{i}((s-1)*step+1,:);
    end
    for i=1:n%sample
        Retriving_Cnt=[s i]
        [v_mem_wta,spikeMat_wta, tVec_wta] = WTA(spikeMat_Coding{i}(:,1:tR),tVec_Coding(:,1:tR),w_stdp_wta);
        spikeMat_wta_plot{i}=spikeMat_wta;
        maxim=find(sum((spikeMat_wta)')==max(sum((spikeMat_wta)')));
        if (maxim==data(i,l+1))
            Crct(s)=Crct(s)+1;
        end
        rate_in{s}(i,:)=(sum(spikeMat_Coding{i}')/tVec_wta(end))*w_stdp_wta';
        rate_out{s}(i,:)=sum(spikeMat_wta')/tVec_wta(end);
    end
end

%% plot Neurons outputs, Correct recognition (Crct), input rate (rate_in) and output rate (rate_out) of neurons
%plot the output firing rate of neurons
figure(1);
for i=1:n % sample
    for j=1:n % neuron
        for s=1:wstep
            rate_O(s,j,i)=rate_out{s}(i,j);
        end
    end
    subplot(2,n/2,i)
    plot(1:wstep,rate_O(:,:,i),'k','LineWidth',2)
    hold on
    plot(1:wstep,rate_O(:,i,i),'-r','LineWidth',2)
    hold on
    xlabel('Time (ms)');ylabel('R.W (MHz)');
    xlim([1 wstep])
end   

%plot the correct pattern recognition versus training time
figure(2)
    plot(1:wstep,Crct/n,'k','LineWidth',2)
    hold on
    xlabel('Time (ms)');ylabel('Original data Correct Recognition Rate');
    xlim([1 wstep]),ylim([0 1.2])
    