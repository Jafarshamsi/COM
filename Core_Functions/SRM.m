function [v_mem,spikeMat_out,tVec_out] = SRM(spikeMat_in,tVec_in,w)
% Spike response model (SRM)
tVec=tVec_in;
nBins = size(tVec_in,2);
nTrials=size(spikeMat_in,1);
v_mem = zeros(1,nBins);
spikeMat_out=zeros(1,nBins);

load('param.mat');

for i=1:nBins
    for j=1:nTrials
        if (spikeMat_in(j,i))
            v_mem=v_mem+w(j)*((exp(-(tVec-tVec(i))/(4*tav))-exp(-(tVec-tVec(i))/tav))).*heaviside(tVec-tVec(i))*amp;
        end
    end
    if (v_mem(i)>0 && v_mem(i)>=vthr && (v_mem(i)-v_mem(i-1)>0))
        v_mem(i)=2;
        v_mem(i+1:end)=0;
        v_mem=v_mem-nu*exp(-(tVec-tVec(i))/4*tav).*heaviside(tVec-tVec(i));
        spikeMat_out(i)=1;
    end
end
tVec_out=tVec_in;
spikeMat_out=logical(spikeMat_out);



