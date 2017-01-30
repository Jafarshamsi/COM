function [w_out,tVec_out] = STDP(spikeMat_in,tVec_in,w_in)
tVec=tVec_in;
nBins = size(tVec_in,2);
nTrials=size(spikeMat_in,1);
v_mem = zeros(1,nBins);
spikeMat_out=zeros(1,nBins);
spiketime_out=zeros(1,nTrials);
spiketime_in=zeros(1,nTrials);
spiketime_in_old=zeros(1,nTrials);

load('param.mat');

npo=0;
npr=0;
npoe=0;
npre=0;
sn=0;
for i=1:nBins
    spiketime_in_old=spiketime_in;
    
    for j=1:nTrials
        if (spikeMat_in(j,i))
            spiketime_in(j)=tVec(i);
            v_mem=v_mem+w_in(end,j)*((exp(-(tVec-tVec(i))/(4*tav))-exp(-(tVec-tVec(i))/tav))).*heaviside(tVec-tVec(i))*amp;
               npo=npo+1;
            if ((spiketime_out(j)>0) &&((tVec(i)-spiketime_out(j))>0) && ((tVec(i)-spiketime_out(j))<3*tav_LTD) )%&& (w_in(end,j)>0)
               npoe=npoe+1;
               w_in(end,j)=w_in(end,j)-A_LTD*exp(-(tVec(i)-spiketime_out(j))/tav_LTD);
%                figure(1)
%                w_d=-A_LTD*exp(-(tVec(i)-spiketime_out(j))/tav_LTD);
%                hold on
%                plot(-(tVec(i)-spiketime_out(j)),w_d,'b.')
               spiketime_out(j)=0;
               if (w_in(end,j)<0)
                   w_in(end,j)=0;
               end
            end
        end
    end
    if (v_mem(i)>0 && v_mem(i)>=vthr && (v_mem(i)-v_mem(i-1)>0))
        v_mem(i)=2;
        v_mem(i+1:end)=0;
        v_mem=v_mem-nu*exp(-(tVec-tVec(i))/tav).*heaviside(tVec-tVec(i));
        spikeMat_out(i)=1;
        spiketime_out(:) = tVec(i);
        sn=sn+1;
        for k=1:nTrials
                 npr=npr+1;
            if ((spiketime_in_old(k)>0) && ((tVec(i)-spiketime_in_old(k))>0) && ((tVec(i)-spiketime_in_old(k))<3*tav_LTP))%&& (w_in(end,k)<wmax) 
                 npre=npre+1;
%                 figure(1)
%                 w_p=A_LTP*exp(-(tVec(i)-spiketime_in_old(k))/tav_LTP);
%                 hold on
%                 plot((tVec(i)-spiketime_in_old(k)),w_p,'r.') 
                w_in(end,k)=w_in(end,k)+A_LTP*exp(-(tVec(i)-spiketime_in_old(k))/tav_LTP);
                 if (w_in(end,k)>wmax)
                   w_in(end,k)=wmax;
                 end
                 spiketime_in_old(k)=0;
                 spiketime_in(k)=0;
            end
            if spiketime_in(k)<tVec(i)
                spiketime_in(k)=0;
            end
        end
    end
    w_in(end+1,:)=w_in(end,:);  
end

% npo
% npr
% npoe
% npre
% sn
w_out=w_in;
tVec_out = tVec_in;