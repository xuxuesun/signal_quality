clc;
clear all;

%preprocessing
sampfreq = 250;
confsig = load('expfile/11.txt');  

confsig1 = smooth(confsig(:,1)-mean(confsig(:,1)));
ppgtt_conf = 1/sampfreq:1/sampfreq:length(confsig1)/sampfreq;

templatesig1 = load('expfile/11template.txt');
templatesig = smooth(templatesig1(:,1)-mean(templatesig1(:,1)));

figure(1);
subplot(2,1,1);
thr1 = -0.75;
thr2 = 0.015;
thr3 = 0.01;
thr4 = 0.01;
[onsetind peakind dicrind secpind mordis] = mydetector11(ppgtt_conf,confsig1,thr1,thr2,thr3,thr4);
h111 = plot(ppgtt_conf(peakind(find(peakind~=-1))), confsig1(peakind(find(peakind~=-1))), 'k^','MarkerSize',10);hold on;
h112 = plot(ppgtt_conf(onsetind), confsig1(onsetind), 'k>','MarkerSize',10);
h113 = plot(ppgtt_conf(dicrind(find(dicrind~=-1))),confsig1(dicrind(find(dicrind~=-1))),'k*','MarkerSize',10);
h114 = plot(ppgtt_conf(secpind(find(secpind~=-1))),confsig1(secpind(find(secpind~=-1))),'ko','MarkerSize',10);
for i=1:length(mordis)
    if mordis(i)==1
       h11=plot(ppgtt_conf(onsetind(i):onsetind(i+1)),confsig1(onsetind(i):onsetind(i+1)),'k','LineWidth',2);
    elseif mordis(i)==-2
       h12=plot(ppgtt_conf(onsetind(i):onsetind(i+1)),confsig1(onsetind(i):onsetind(i+1)),'k:','LineWidth',3);
    end
end
legend([h111 h112 h113 h114 h11 h12],'Peak','Onset','Dicrotic Notch','Dicr Peak','Level1: good','Level4: corrupted');
xlabel('Time(seconds)');
ylabel('Amplitude');
title('(a) stage1 analysis result');

mymorflag = mordis;
for i=1:length(onsetind)-1
    if mymorflag(i)>=0
        [myorigcost(i) mydercost(i)]= dtwcost(confsig1(onsetind(i):onsetind(i+1)),templatesig);
        if  mydercost(i) > 1 & mymorflag(i)==1     %origcost(i) > 2 ||
           mymorflag(i)=0;
        end
        if mymorflag(i)==0
           % plot(ppgtt_conf(onsetind(i):onsetind(i+1)),confsig1(onsetind(i):onsetind(i+1)),'k');
        end
    end
end

initf = 'expfile/11init.txt';%rnd51init.txt,rnd21init.txt
mysigqly = kalesti1(ppgtt_conf,initf,confsig1,peakind,onsetind,dicrind,secpind,222,mymorflag);
subplot(2,1,2);
plot(ppgtt_conf(1:onsetind(1)),confsig1(1:onsetind(1)),'k','LineWidth',2);hold on;
plot(ppgtt_conf(onsetind(end):length(confsig1)),confsig1(onsetind(end):length(confsig1)),'k','LineWidth',2);
for i=1:length(onsetind)-1
   if mysigqly(i)==-2
       h4 = plot(ppgtt_conf(onsetind(i):onsetind(i+1)),confsig1(onsetind(i):onsetind(i+1)),'k:','LineWidth',3);
   elseif mysigqly(i)==-1
       h3 = plot(ppgtt_conf(onsetind(i):onsetind(i+1)),confsig1(onsetind(i):onsetind(i+1)),'k-.','LineWidth',3);
   elseif mysigqly(i)==0
       h2 = plot(ppgtt_conf(onsetind(i):onsetind(i+1)),confsig1(onsetind(i):onsetind(i+1)),'k--','LineWidth',3);
   else
       h1 = plot(ppgtt_conf(onsetind(i):onsetind(i+1)),confsig1(onsetind(i):onsetind(i+1)),'k','LineWidth',2);
   end
end
h5 = plot(ppgtt_conf(onsetind), confsig1(onsetind), 'k>','MarkerSize',10);
legend([h4 h3 h2 h1 h5],'Level4: total corrupted','Level3: invalid features','Level2: valid features, bad morphology','Level1: good','onset');
xlabel('Time(seconds)');
ylabel('Amplitude');
title('(b) final assessment result');
hold off;
grid on;
axis auto fill;