% This simulation is itrated corresponding to table 1 and the output are saved in files (.Mat).
% This code inclues :
%                   1. Generation of the Dataset;
%                   2. Adding noise to the Dataset;
%                   3. Coding the orginal dataset based on Poisson distribution;
%                   4. Coding the noisy data based on Poisson distribution;
%                   5. STDP training of a COM;
%                   6. Processing of the COM for orginal data;
%                   7. Processing of the COM for noisy data;
%                   8. Ploting the output of neurons for COM vs WTAs;
%                   9. Ploting the Recognition rate of orginal and noisy.
% How to use:
%                   1. set the parameters 
%                   2. Set ERS (0=fig5a, 0.25=fig5b, 0.5=fig5c)
%                   3. RUN
%% CLEAR

clear
clc

%% parameters
for lnNM_num=1:10
%lnNM is corresponding to Table 1 in paper
lnNM=[8 8 4 10;10 10 4 10;12 12 4 10;14 14 4 10;16 8 4 10;16 16 4 10;16 32 4 10;32 16 4 10;32 32 4 10;32 64 4 10]; 
nn=8; % Number of noise (0.1 to 0.8)
tT=20;% Training time
tR=300;% Retrival time
alpha=0.4;% alpha is related to wext in Eq. 3
ERS = 0;% Erasing the message (0, 0.25, 0.5)

mu_weight_init=0.2;% Mu initial weight
sigma_weight_init=0.01;% Sigma initial weight

%% simulation

l=lnNM(lnNM_num,1); % length of pattern
n=lnNM(lnNM_num,2); % Number of patterns = Number of neuron
N=lnNM(lnNM_num,3); % Number of patterns in a message = Number of WTAs
M=lnNM(lnNM_num,4);% Number of messages


%% Dataset;OD;MNIST,RAND
        sl=sqrt(l);
        for j=1:N
            for i=1:n
                data{j}(i,1:l)=(0.8+0.1*rand(1,l)).*(-floor(-full(sprand(1,l,0.5))));% generate random dataset through Gaussian distribution
                data{j}(i,l+1)=i;
            end
        end

%Generation of Messages
    for m=1:M
        message(m,:)=randperm(n,N);
    end
%% Adding noise to Dataset;OD;MNIST,RAND

 gn=[0.1:0.1:0.1*nn];
 for noise=1:nn
     for i=1:N
        gnoise=gn(noise)*(2*rand(n,l)-1);
        ndata{noise,i}(:,1:l) =data{i}(:,1:l)+gnoise;
        ndata{noise,i}(:,l+1) =data{i}(:,l+1);
        ndata{noise,i}(find( ndata{noise,i}(:,1:l)<0))=0;
        ndata{noise,i}(find( ndata{noise,i}(:,1:l)>1))=1;
     end
 end


%% Coding 

for j=1:N
    for i=1:n
        Coding_Cnt=[j i]
        [spikeMat_Coding1, tVec_Coding] = Coding(tT,data{j}(i,1:l));
        spikeMat_Coding{j}{i}(:,:)=spikeMat_Coding1;
    end
end

 %% Coding Noisy data

for noise=1:nn
    for j=1:N
        for i=1:n
            nCoding_Cnt=[noise i]
            [spikeMat_Coding1, tVec_Coding] = Coding(tT,ndata{noise,j}(i,1:l));
            nspikeMat_Coding{noise}{j}{i}(:,:)=spikeMat_Coding1;
        end
    end
end

%% Training

w_init_STDP=mu_weight_init + sigma_weight_init.*randn(1, l);

for j=1:N
    for i=1:n
        Training_Cnt=[j,i]
        [w_stdp1,tVec_stdp] = STDP(spikeMat_Coding{j}{i},tVec_Coding,w_init_STDP);
        w_stdp{j}(i,:)=w_stdp1(end,:);
    end
end

for m=1:M
    for N1=1:N
        for N2=1:N
            if (N1 ~= N2)
                w_hebb(N1,message(m,N1),N2,message(m,N2)) = 1/(alpha*N);
            end
        end
    end
end

%% Earasing data
EN=ERS*N; % Erasing the message (0, 0.25, 0.5)
E1=randperm(N,2);
for k=1:EN
    for e=1:n
        spikeMat_Coding{E1(k)}{e}=spikeMat_Coding{E1(k)}{e}*0;
    end
    for en=1:nn
        for enn=1:n
            nspikeMat_Coding{en}{E1(k)}{enn}=nspikeMat_Coding{en}{E1(k)}{enn}*0;
        end
    end
end
%% Processing Orginal data: Neurons outputs COM & individual WTAs, Correct recognition (CrctCOM & CrctWTA)
CrctCOM=0;
CrctWTA=0;
for j=1:M%select message
        Retriving_Cnt=j
        for k=1:N
            spikeMat_Coding_M1=spikeMat_Coding{k}(message(j,k));%message spike matrix
            spikeMat_Coding_M(k) = {spikeMat_Coding_M1{1}(:,1:tR)};
        end
        [v_mem_WTA,spikeMat_WTA, tVec_WTA] = COM(spikeMat_Coding_M,tVec_Coding(:,1:tR),w_stdp,w_hebb*0);
        spikeMat_WTAs_plot{j}=spikeMat_WTA;
        [v_mem_COM,spikeMat_COM, tVec_COM] = COM(spikeMat_Coding_M,tVec_Coding(:,1:tR),w_stdp,w_hebb);
        spikeMat_COM_plot{j}=spikeMat_COM;
        for i=1:N
            maxim_WTA=find(sum((spikeMat_WTA{i})')==max(sum((spikeMat_WTA{i})')))
            if (maxim_WTA==message(j,i))
                CrctWTA=CrctWTA+1;
            end
            maxim_COM=find(sum((spikeMat_COM{i})')==max(sum((spikeMat_COM{i})')));
            if (maxim_COM==message(j,i))
                CrctCOM=CrctCOM+1;
            end
        end
        CrctCOM
        CrctWTA
end

%% Processing Noisy data: Correct recognition (Crct)

nCrctCOM=zeros(1,nn);
nCrctWTA=zeros(1,nn);
for noise=1:nn
    for j=1:M%select message
        Retriving_Cnt=[noise j]
        for k=1:N
            spikeMat_Coding_M1=nspikeMat_Coding{noise}{k}(message(j,k));%message spike matrix
            spikeMat_Coding_M(k) = {spikeMat_Coding_M1{1}(:,1:tR)};
        end
        [v_mem_COM,spikeMat_COM, tVec_COM] = COM(spikeMat_Coding_M,tVec_Coding(:,1:tR),w_stdp,w_hebb);
        [v_mem_WTA,spikeMat_WTA, tVec_WTA] = COM(spikeMat_Coding_M,tVec_Coding(:,1:tR),w_stdp,w_hebb*0);
        for i=1:N
            maxim=find(sum((spikeMat_COM{i})')==max(sum((spikeMat_COM{i})')));
            if (maxim==message(j,i))
                nCrctCOM(noise)=nCrctCOM(noise)+1;
            end
        end
        for i=1:N
            maxim=find(sum((spikeMat_WTA{i})')==max(sum((spikeMat_WTA{i})')));
            if (maxim==message(j,i))
                nCrctWTA(noise)=nCrctWTA(noise)+1;
            end
        end
        nCrctCOM
        nCrctWTA
    end
end
save(strcat('COMvsWTA',num2str(ERS),num2str(lnNM_num),'.mat'));


%% Clear
clear % clear everything for next simulation

end
%% plot Neurons outputs, Correct recognition (Crct)
% This part of the code is used to plot the curves

CrctCOMFinal=0;
CrctWTAFinal=0;
TN=10;
for i=1:TN
load(strcat('COMvsWTA0',num2str(i),'.mat'));%COMvsWTA0=fig5a,COMvsWTA25=fig5b, or COMvsWTA50=fig5b
Crct_COM = [CrctCOM nCrctCOM];
Crct_WTA = [CrctWTA nCrctWTA];
CrctCOMFinal=Crct_COM+CrctCOMFinal;
CrctWTAFinal=Crct_WTA+CrctWTAFinal;
end
CrctCOMFinal=CrctCOMFinal/(40);
CrctWTAFinal=CrctWTAFinal/(40);

figure(1)
    plot(0:0.1:0.8,CrctCOMFinal/TN,'r','LineWidth',3)
    hold on
    plot(0:0.1:0.8,CrctWTAFinal/TN,'--b','LineWidth',3)
    xlabel('Noise(\sigma)');ylabel('Correct Recognition Rate');
    xlim([0 0.8]),ylim([0 1.2])