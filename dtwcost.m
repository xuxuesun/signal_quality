function [ origco derco ] = dtwcost( sig1,templsig )
%DTWCOST Summary of this function goes here
%   Detailed explanation goes here
%[seq1,seq2,distmatr] = dp(ddtw_dist_calc(sig1,templsig));
[seq1,seq2,distmatr] = dp(dtw_dist_calc(sig1,templsig));
[dseq1,dseq2,ddistmatr] = dp(ddtw_dist_calc(sig1,templsig));
%calculate the cost!!
origco = 0;
derco = 0;
for i=1:length(seq1)
    origco = origco + distmatr(seq1(i),seq2(i));
end
for i=1:length(dseq1)
    derco = derco + ddistmatr(dseq1(i),dseq2(i));
end
end

function [ seqQ,seqC,cost ] = classicdtw( sigq, sigc )
DTWval = dtw_dist_calc(sigq,sigc);


%subplot(4,1,2);
%sigq= load('template11.txt');
%plot();
end




%template:1000Hz

function mat = wddtw_dist_calc(serq,serc)
weival= zeros(length(serq),length(serc));
wmax = 50;
%weival = wmax/(1+e^(-0.25*(median(serq)));
end

function mat = ddtw_dist_calc(serq,serc)
mat = zeros(length(serq),length(serc));
estiq = zeros(length(serq));
estic = zeros(length(serc));
mat = zeros(length(serq),length(serc));
for i=2:length(serq)-1
       estiq(i) = (serq(i)-serq(i-1)+(serq(i+1)-serq(i-1))/2)/2;
end
estiq(1) = estiq(2);estiq(length(serq)) =estiq(length(serq)-1);
for i=2:length(serc)-1
       estic(i) = (serc(i)-serc(i-1)+(serc(i+1)-serc(i-1))/2)/2;
end
estic(1) = estic(2);estic(length(serc)) =estic(length(serc)-1);
for i= 1:length(estiq)
    for j=1:length(estic)
        mat(i,j) = (estiq(i)-estic(j))^2;
    end
end
end

function mat = dtw_dist_calc(serq,serc)
mat = zeros(length(serq),length(serc));
for i=1:length(serq)
    for j=1:length(serc)
        mat(i,j)=(serq(i)-serc(j))^2;
    %    if i==1 && j~=1
    %        mat(i,j) = mat(i,j-1)+mat(i,j);
    %    elseif i~=1 && j==1
    %        mat(i,j) = mat(i-1,j)+mat(i,j);  
    %    elseif i~=1&&j~=1
    %        mat(i,j) = min(mat(i,j-1),mat(i-1,j),mat(i-1,j-1))+mat(i,j);    
    %    end
    end
end
end
