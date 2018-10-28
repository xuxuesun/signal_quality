function [sigqual] = kalmestimate(timespan,initfile,estisig,peakpos,onsetpos,dicrpos,secppos,figname,morass)
initsig1 = load(initfile);
initsig = smooth(initsig1(:,1)-mean(initsig1(:,1)));
thr1 = -0.75;
thr2 = 0.015;
thr3 = 0.015;
thr4 = 0.003;
timeinit = 1/250:1/250:length(initsig)/250;
[onsetposinit,peakposinit,dicroninit]=mydelineator(initsig,250);

figure(figname);
subplot(2,1,1);
plot(timeinit, initsig);hold on;
plot(timeinit(onsetposinit), initsig(onsetposinit), 'm>');
plot(timeinit(peakposinit(find(peakposinit~=-1))), initsig(peakposinit(find(peakposinit~=-1))), 'r^');
%plot(timeinit(dicroninit(find(dicroninit~=-1))),initsig(dicroninit(find(dicroninit~=-1))),'g*');
%plot(timeinit(secpinit(find(secpinit~=-1))),initsig(secpinit(find(secpinit~=-1))),'k.');
hold off;
subplot(2,1,2);
plot(timespan, estisig);hold on;
plot(timespan(onsetpos), estisig(onsetpos), 'm>');
plot(timespan(peakpos(find(peakpos~=-1))), estisig(peakpos(find(peakpos~=-1))), 'r^');
plot(timespan(dicrpos(find(dicrpos~=-1))), estisig(dicrpos(find(dicrpos~=-1))), 'g*');
if isempty(secppos)
    plot(timespan(secppos(find(secppos~=-1))), estisig(secppos(find(secppos~=-1))), 'k.');
end
%initiate kalman parameter and observation
meas(1,1:length(onsetpos)-1) = diff(timespan(onsetpos(1:length(onsetpos))));        %onset duration
for i=1:length(onsetpos)-2
   if peakpos(i+1) ~= -1 & peakpos(i)~=-1
       meas(2,i) = diff(timespan(peakpos(i:i+1)));       %pp interval
   else
       meas(2,i) = -1;
   end
   if peakpos(i)~=-1
       meas(3,i) = estisig(peakpos(i))-estisig(onsetpos(i));   
   else
       meas(3,i)=-1;
   end
end
meas(2,length(onsetpos)-1) = 0.678;
if peakpos(length(onsetpos)-1)~=-1
    meas(3,length(onsetpos)-1) = estisig(peakpos(length(onsetpos)-1))-estisig(onsetpos(length(onsetpos)-1));
else
    meas(3,length(onsetpos)-1) =-1;
end
meas(4,1:length(onsetpos)-1)=estisig(onsetpos(1:length(onsetpos)-1));   %onset y axis location

initval(1) = mean(diff(timespan(onsetposinit)));
initval(2) = mean(diff(timespan(peakposinit)));
mininitlen = min(length(onsetposinit),length(peakposinit));
initval(3) = mean(initsig(peakposinit(1:mininitlen))-initsig(onsetposinit(1:mininitlen)));
initval(4)=mean(initsig(onsetposinit(1:mininitlen)));
covr1 = cov(diff(timespan(onsetposinit)));
covr2 = cov(diff(timespan(peakposinit)));
covr3 = cov(initsig(peakposinit(1:mininitlen))-initsig(onsetposinit(1:mininitlen)));
covr4 = cov(initsig(onsetposinit(1:mininitlen)));
kalmesti.R=[covr1 0 0 0
    0 covr2 0 0
    0 0 covr3 0
    0 0 0 covr4];
kalmesti.A = eye(4);
kalmesti.x = initval';
kalmesti.H=eye(4);
kalmesti.P = inv(kalmesti.H)*kalmesti.R*inv(kalmesti.H');
kalmesti.u=0;
kalmesti.B=0;
kalmesti.Q=zeros(4);
%kalmesti.Q=eye(3);
%estimation results
ppinterv_thred1 = 1;
ppinterv_thred2 = 0.4;
segdura_thred = 0.85;
piresi_thred = initval(2)/4;
odresi_thred = initval(1)/4;
heresi_thred = initval(3)/2;
ohresi_thred = 0.4;
sigqual = morass;    %peak-to-peak interval abnormal/normal:duration abnormal/normal:height abnormal/normal      1/0:1/0:1/0
piass = zeros(size(meas(1,:)));
peakass = zeros(size(meas(1,:))+1);
for m=1:length(meas)
    if morass(m) >=0            %segment need to be estimated
       kalmesti(end).z=meas(:,m);
       [kalmesti(end+1) xpred xresi]=kalmanf(kalmesti(end),piresi_thred,odresi_thred,heresi_thred,ohresi_thred);
       
       if meas(2,m)==-1
          piass(m)= -1;
          if peakpos(m) ~= -1
              peakass(m)=1;
          else
              peakass(m)=-1;
          end
      %   kalmesti(end).x(2,:) = xpred(2,:);
      %   tempP = kalmesti(end-1).A * kalmesti(end-1).P * kalmesti(end-1).A' + kalmesti(end-1).Q;
      %   kalmesti(end).P(2,:) = tempP(2,:);
       else
           if m==1
               if abs(xresi(2,:))>piresi_thred && (morass(m)+morass(m+1)<2)
                   peakass(m) = 0;
                   piass(m) = 0;
                   plot(timespan(peakpos(m)),estisig(peakpos(m)),'k.');
               else
                   peakass(m) = 1;
                   piass(m) = 1;
               end               
           else
               if abs(xresi(2,:))>piresi_thred
                   if m~=length(meas)
                       if morass(m)+morass(m+1)<2
                           piass(m)=0;
                           if piass(m-1) ==0
                               peakass(m)=0;
                               plot(timespan(peakpos(m)),estisig(peakpos(m)),'k.');
                           end
                       else
                           piass(m)=1;
                           if piass(m-1) ==1
                               peakass(m)=1;
                           else
                               bar(timespan(onsetpos(m)),max(estisig),0.05);
                           end
                       end
                   else
                       if morass(m)<1
                           peakass(m) = 0;
                           if peakpos(m)~=-1
                               plot(timespan(peakpos(m)),estisig(peakpos(m)),'k.');
                           end
                       else
                           peakass(m) =1;
                           piass(m)=1;
                       end
                   end
               else
                   piass(m) = 1;
                   if piass(m-1)==1
                       peakass(m) =1;
                   else
                       %label new start for calculating peak interval 
                       bar(timespan(onsetpos(m)),max(estisig),0.05);
                       peakass(m) =1;
                       %arrow3([timespan(onsetpos(m)) max(estisig)],[timespan(onsetpos(m))+0.1 max(estisig)],'e');
                   end
                   if m == length(meas)-1
                       peakass(m)=1;
                       peakass(m+1) = 1;
                   end
               end
           end
       end
     
       if abs(xresi(1,:))>odresi_thred & morass(m)<1
           sigqual(m)=-1;
       else
           if (abs(xresi(3,:))> heresi_thred & meas(3,m)~=-1) || abs(xresi(4,:))>ohresi_thred
               if m+1<=length(peakpos) && peakpos(m+1)~=-1
                   plot(timespan(onsetpos(m):peakpos(m)), estisig(onsetpos(m):peakpos(m)),'g--');
               elseif m==length(peakpos)
                   plot(timespan(onsetpos(m):length(estisig)), estisig(onsetpos(m):length(estisig)),'g--');
               end
               sigqual(m)=0;
           %else
           %    sigqual(m)=1;
           end
       end
       
       if piass(m)==0 | sigqual(m)<1
           if sigqual(m)<0
               plot(timespan(onsetpos(m):onsetpos(m+1)), estisig(onsetpos(m):onsetpos(m+1)),'r.-');
           else
               plot(timespan(onsetpos(m):onsetpos(m+1)), estisig(onsetpos(m):onsetpos(m+1)),'y.-');               
           end
         %kalmesti(end) = kalmesti(end-1);
         %kalmesti(end).x = xpred;
         %kalmesti(end).P = kalmesti(end).A * kalmesti(end).P * kalmesti(end).A' + kalmesti(end).Q;
       else
           %plot(timespan(onsetpos(m))+kalmesti(end).x(1,:),kalmesti(end).x(4,:),'k>');
           stem(timespan(onsetpos(m))+kalmesti(end).x(1,:),-2);
           %stem(timespan(onsetpos(m))+xpred(1,:),-2);
           if m ~=length(meas) && piass(m)==1
               %plot(timespan(peakpos(m))+kalmesti(end).x(2,:),estisig(onsetpos(m+1))+kalmesti(end).x(3,:),'k^');
               plot(timespan(peakpos(m))+xpred(2,:),estisig(onsetpos(m+1))+xpred(3,:),'k^');
           end
       end
    else
        plot(timespan(onsetpos(m):onsetpos(m+1)), estisig(onsetpos(m):onsetpos(m+1)),'r');
    end
end
pqind = find(peakass ==0);
flagind = find(sigqual>=0);
sigqual(intersect(pqind,flagind)) = -1;
hold off;
grid on
axis auto fill;
end