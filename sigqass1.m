clc;
clear all;

%preprocessing
sampfreq = 250;
confsig = load('expfile/11.txt');        %rnd51.txt,rnd21.txt

confsig1 = smooth(confsig(:,1)-mean(confsig(:,1)));
ppgtt_conf = 1/sampfreq:1/sampfreq:length(confsig1)/sampfreq;

templatesig1 = load('expfile/11template.txt');       %rnd51template.txt,rnd21template.txt
templatesig = smooth(templatesig1(:,1)-mean(templatesig1(:,1)));

figure(1);
subplot(2,1,1);
%segmentation & feature points detection
[onsetp,peakp,dicron]=mydelineator(confsig1,250);
plot(ppgtt_conf, confsig1);hold on;
plot(ppgtt_conf(onsetp), confsig1(onsetp), 'm>');
plot(ppgtt_conf(peakp), confsig1(peakp), 'k^');
plot(ppgtt_conf(dicron), confsig1(dicron), 'g*');
%legend('PPG Waveforms', 'Onset','Peak', 'Dicron');
xlabel('Time');
ylabel('Signal & Delineation'); 
grid on;
axis auto fill;
%Template dynamic time warping calculation
morflag = [];
for i=1:length(onsetp)-1
    [origcost(i) dercost(i)]= dtwcost(confsig1(onsetp(i):onsetp(i+1)),templatesig);
    if dercost(i) > 1
        if origcost(i)>200
            morflag = [morflag -2];
            plot(ppgtt_conf(onsetp(i):onsetp(i+1)),confsig1(onsetp(i):onsetp(i+1)),'r');
        elseif origcost(i)>100
            morflag = [morflag -1];
            plot(ppgtt_conf(onsetp(i):onsetp(i+1)),confsig1(onsetp(i):onsetp(i+1)),'k');
        else
            morflag = [morflag 0];
        end
    else
        morflag = [morflag 1];
    end        
end
hold off;

subplot(2,1,2);
thr1 = -0.75;
thr2 = 0.015;
thr3 = 0.01;
thr4 = 0.01;
[onsetind peakind dicrind secpind mordis] = mydetector11(ppgtt_conf,confsig1,thr1,thr2,thr3,thr4);
plot(ppgtt_conf, confsig1);hold on;
plot(ppgtt_conf(peakind(find(peakind~=-1))), confsig1(peakind(find(peakind~=-1))), 'k^');
plot(ppgtt_conf(onsetind), confsig1(onsetind), 'm>');
plot(ppgtt_conf(dicrind(find(dicrind~=-1))),confsig1(dicrind(find(dicrind~=-1))),'g*');
plot(ppgtt_conf(secpind(find(secpind~=-1))),confsig1(secpind(find(secpind~=-1))),'r.');
for i=1:length(mordis)
    if mordis(i)<0
        plot(ppgtt_conf(onsetind(i):onsetind(i+1)),confsig1(onsetind(i):onsetind(i+1)),'r');
    end
end
%Template dynamic time warping calculation
mymorflag = mordis;
for i=1:length(onsetind)-1
    if mymorflag(i)>=0
        [myorigcost(i) mydercost(i)]= dtwcost(confsig1(onsetind(i):onsetind(i+1)),templatesig);
        if  mydercost(i) > 1 & mymorflag(i)==1     %origcost(i) > 2 ||
           mymorflag(i)=0;
        end
        if mymorflag(i)==0
            plot(ppgtt_conf(onsetind(i):onsetind(i+1)),confsig1(onsetind(i):onsetind(i+1)),'k');
        end
    end
end
grid on;
axis auto fill;
hold off;

%%kalman prediction&estimation
%quasi-periodic duration and peak interval estimation
initf = 'expfile/11init.txt';%rnd51init.txt,rnd21init.txt
sigqly = kalmestimate(ppgtt_conf,initf,confsig1,peakp,onsetp,dicron,[],111,morflag);
mysigqly = kalmestimate(ppgtt_conf,initf,confsig1,peakind,onsetind,dicrind,secpind,222,mymorflag);
figure(611);
subplot(2,1,1);
plot(ppgtt_conf, confsig1,'k');hold on;
for i=1:length(onsetp)-1
   if sigqly(i)==-2
       plot(ppgtt_conf(onsetp(i):onsetp(i+1)),confsig1(onsetp(i):onsetp(i+1)),'r');
   elseif sigqly(i)==-1
       plot(ppgtt_conf(onsetp(i):onsetp(i+1)),confsig1(onsetp(i):onsetp(i+1)),'r.');
   elseif sigqly(i)==0
       plot(ppgtt_conf(onsetp(i):onsetp(i+1)),confsig1(onsetp(i):onsetp(i+1)),'y.');
   else
       plot(ppgtt_conf(onsetp(i):onsetp(i+1)),confsig1(onsetp(i):onsetp(i+1)),'b.');
   end
end
hold off;

figure(612);
subplot(2,1,2);
plot(ppgtt_conf(1:onsetind(1)),confsig1(1:onsetind(1)),'k');hold on;
plot(ppgtt_conf(onsetind(end):length(confsig1)),confsig1(onsetind(end):length(confsig1)),'k');
for i=1:length(onsetind)-1
   if mysigqly(i)==-2
       h4 = plot(ppgtt_conf(onsetind(i):onsetind(i+1)),confsig1(onsetind(i):onsetind(i+1)),'r');
   elseif mysigqly(i)==-1
       h3 = plot(ppgtt_conf(onsetind(i):onsetind(i+1)),confsig1(onsetind(i):onsetind(i+1)),'r--');
   elseif mysigqly(i)==0
       h2 = plot(ppgtt_conf(onsetind(i):onsetind(i+1)),confsig1(onsetind(i):onsetind(i+1)),'y.');
   else
       h1 = plot(ppgtt_conf(onsetind(i):onsetind(i+1)),confsig1(onsetind(i):onsetind(i+1)),'b');
   end
end
h5 = plot(ppgtt_conf(onsetind), confsig1(onsetind), 'm>');
legend([h4 h3 h2 h1 h5],'Level4: total corrupted','Level3: invalid features','Level2: valid features, bad morphology','Level1: good','onset');
xlabel('Time(seconds)');
ylabel('Amplitude');
hold off;
grid on;
axis auto fill;