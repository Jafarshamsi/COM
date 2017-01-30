function [v_mem,spikeMat_out,tVec_out] = COM(spikeMat_in,tVec_in,wI2C,wC2C)
tVec=tVec_in;
nBins = size(tVec_in,2);
N=size(spikeMat_in,2);% # of minicolumn
nTrials=size(spikeMat_in{1},1);% # of trial in each minicolumn
nNeuron=size(wI2C{1},1);% # of neuron in each mincolumn

v_mem =cell(N,1);v_mem(:,1)={zeros(nNeuron,nBins)};
spikeMat_out=cell(N,1);spikeMat_out(:,1)={zeros(nNeuron,nBins)};

load('param.mat');

spikeMat_out=cell(N,1);
for ii=1:N
    spikeMat_out{ii}=zeros(nNeuron,nBins);
end
for timer=1:nBins-1
    for m=1:N% v_mem of neuron from own column
        for k=1:nNeuron
            for j=1:nTrials
                if (spikeMat_in{m}(j,timer))
                    v_mem{m}(k,:)=v_mem{m}(k,:)+wI2C{m}(k,j)*((exp(-(tVec-tVec(timer))/(4*tav))-exp(-(tVec-tVec(timer))/tav))).*heaviside(tVec-tVec(timer))*amp;
                end
            end
        end
    end    

    for m=1:N% spike out of neurons and reseting the neuron
        for k=1:nNeuron 
            if (v_mem{m}(k,timer)>0 && v_mem{m}(k,timer)>=vthr && (v_mem{m}(k,timer)-v_mem{m}(k,timer-1)>0))
                v_mem{m}(k,timer)=2;
                v_mem{m}(k,timer+1:end)=0;
                %v_mem{m}(k,:)=v_mem{m}(k,:)-nu*exp(-(tVec-tVec(timer))/tav).*heaviside(tVec-tVec(timer));
                spikeMat_out{m}(k,timer)=1;
            end
        end
    end
    for m=1:N
        for k=1:nNeuron
            if (spikeMat_out{m}(k,timer)==1)
                v_mem{m}([1:k-1 k+1:end],timer+1:end)=0;
                %v_mem{m}([1:k-1 k+1:end],timer+1:end)=v_mem{m}([1:k-1 k+1:end],:)+repmat((-2)*((exp(-(tVec-tVec(i+1))/(4*tav))-exp(-(tVec-tVec(i+1))/(tav)))).*heaviside(tVec-tVec(i+1))*amp,nNeuron-1,1);
            end
        end
    end
    for m=1:N% v_mem of neuron from other column
        for k=1:nNeuron
            for mm=1:N
                 for kk=1:nNeuron 
                    if (spikeMat_out{m}(k,timer))
                        v_mem{mm}(kk,:)=v_mem{mm}(kk,:)+wC2C(m,k,mm,kk)*((exp(-(tVec-tVec(timer+1))/(4*tav))-exp(-(tVec-tVec(timer+1))/(tav)))).*heaviside(tVec-tVec(timer+1))*amp;
                    end
                 end
            end
        end
    end
end
tVec_out=tVec_in;
for m=1:N
    spikeMat_out{m}=logical(spikeMat_out{m});
end



