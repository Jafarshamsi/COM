function [v_mem,spikeMat_out,tVec_out] = WTA(spikeMat_in,tVec_in,w)
tVec=tVec_in;
nBins = size(tVec_in,2);
nTrials=size(spikeMat_in,1);
nNeuron=size(w,1);
v_mem = zeros(nNeuron,nBins);
spikeMat_out=zeros(nNeuron,nBins);

load('param.mat');

for i=1:nBins
    for k=1:nNeuron
        for j=1:nTrials
            if (spikeMat_in(j,i))
                v_mem(k,:)=v_mem(k,:)+w(k,j)*((exp(-(tVec-tVec(i))/(4*tav))-exp(-(tVec-tVec(i))/tav))).*heaviside(tVec-tVec(i))*amp;
            end
        end
    end
    for k=1:nNeuron 
        if (v_mem(k,i)>0 && v_mem(k,i)>=vthr && (v_mem(k,i)-v_mem(k,i-1)>0))
            v_mem(k,i)=2;
            v_mem(k,i+1:end)=0;
            %v_mem(k,:)=v_mem(k,:)-nu*exp(-(tVec-tVec(i))/tav).*heaviside(tVec-tVec(i));
            spikeMat_out(k,i)=1;
        end
    end
    for k=1:nNeuron
        if (spikeMat_out(k,i)==1 && i<size(tVec,2))
            v_mem([1:k-1 k+1:end],i+1:end)=0;
            %v_mem([1:k-1 k+1:end],:)=v_mem([1:k-1 k+1:end],:)+repmat(-nu*exp(-(tVec-tVec(i+1))/(4*tav)).*heaviside(tVec-tVec(i+1)),nNeuron-1,1);
        end
    end
end
tVec_out=tVec_in;
spikeMat_out=logical(spikeMat_out);



