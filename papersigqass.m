clc;
clear all;

%preprocessing
sampfreq = 250;
confsig = load('expfile/rnd51.txt');        %11,rnd21,21,51,rnd51

confsig1 = smooth(confsig(:,1)-mean(confsig(:,1)));
ppgtt_conf = 1/sampfreq:1/sampfreq:length(confsig1)/sampfreq;

templatesig1 = load('expfile/rnd51template.txt');   
templatesig = smooth(templatesig1(:,1)-mean(templatesig1(:,1)));

figure(1);
subplot(2,1,1);
%segmentation & feature points detection
[onsetp,peakp,dicron]=mydelineator(confsig1,500);
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

initf = 'expfile/rnd51init.txt';%rnd51init.txt,rnd21init.txt
%morflag(:)=1;
sigqly = kalmestimate(ppgtt_conf,initf,confsig1,peakp,onsetp,dicron,[],111,morflag);
figure(611);
subplot(2,1,1);
plot(ppgtt_conf, confsig1,'k');hold on;
plot(ppgtt_conf(onsetp),confsig1(onsetp),'k>');
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
grid on;
axis auto fill;