function [ onsind peaind dicrind secpind freqdis] = mydetector11( timespan,detsig,ampthr,indthr,peakthr,dicrthr)
medwinsiz1 = 125;
ressig = detsig-medfilt1(detsig,medwinsiz1);

ind = intersect(find(ressig==0),union(find(diff(ressig)==0),find(diff(ressig)==0)+1));

ind2 = union(find(diff(ind)~=1)+1,find(diff(ind)~=1));
ind2 = ind(ind2);
indleft = [];
indright = [];
for i=1:length(ind2)
   tmpind = ind2(i);
   if ressig(tmpind-1)==0 && ressig(tmpind+1)~=0
       indleft = [indleft tmpind];
   elseif ressig(tmpind-1)~=0 && ressig(tmpind+1)==0
       indright = [indright tmpind];
   end
end

onsind = [];
if indleft(1)>150
    onsind = [1 onsind];
else
    firstdur = 1:indleft(1);
    thredind3 = find(abs(ressig(firstdur))<=indthr);
    ressig(firstdur(thredind3)) = 0;
    diffsig3 = diff(ressig(firstdur));
    firstind = find(diffsig3(1:length(diffsig3)-1).*diffsig3(2:length(diffsig3))<0 & diffsig3(1:length(diffsig3)-1)<0);
    if ~isempty(firstind)
        onsind = [firstdur(firstind+1) onsind];
    end
end
for i=1:length(indright)
    indrang = indleft(i)+1:indright(i)-1;
    ampthrind = find(ressig(indrang)<=ampthr);   
    tmpind = find(ressig(indrang)==0);
    candind = -1;
    if isempty(tmpind) && ressig(indleft(i)+1)-ressig(indleft(i))<0 && ressig(indright(i))-ressig(indright(i)-1)>0  %|| ~isempty(ampthrind)
        [minval minind] = min(ressig(indrang));
        candind = indrang(minind);        
    else
            thredind = find(abs(ressig(indrang))<=indthr);
            ressig(indrang(thredind)) = 0;
            newleftind =  intersect(find(ressig(indrang)==0),find(diff(ressig(indrang)~=0)));
            if length(newleftind)==0
                newleftind = indleft(i);
            else
                newleftind = indrang(newleftind);
            end
            newrightind = intersect(find(ressig(indrang)==0),find(diff(ressig(indrang)~=0))+1);
            if length(newrightind)==0
                newrightind = indright(i);
            else
                newrightind = indrang(newrightind);
            end
            if length(newleftind)==1 && length(newrightind)==1
                tmpind2 = find(ressig(newleftind+1:newrightind-1)==0);
                 if isempty(tmpind2) && ressig(newleftind+1)-ressig(newleftind)<0 && ressig(newrightind)-ressig(newrightind-1)>0 && ~isempty(thredind)
                    [minval minind] = min(ressig(newleftind+1:newrightind-1));
                    newrang = newleftind+1:newrightind-1;
                    candind = newrang(minind);
                end
            end
    end    
    if candind ~= -1
        if detsig(candind)-detsig(candind-1)<0 && detsig(candind+1)-detsig(candind)>0 %|| ~isempty(ampthrind)
            onsind = [onsind candind];
        else
            %search in 5samples window
            serwin1 = candind-3:candind+2;
            diffsig1 = diff(detsig(serwin1));
            newcanind = find(diffsig1(1:4).*diffsig1(2:5)<0 & diffsig1(1:4)<0);
            if ~isempty(newcanind)
                onsind = [onsind serwin1(newcanind+1)];
            end
        end
    end
end

if indright(end) < length(ressig)
    lastdur = indright(length(indright)):length(ressig);
    thredind2 = find(abs(ressig(lastdur))<=indthr);
    ressig(lastdur(thredind2)) = 0;
    lastseg = ressig(lastdur);
    diffsig2 =  diff(lastseg);
    lastind = find(diffsig2(1:length(diffsig2)-1).*diffsig2(2:length(diffsig2))<0 & diffsig2(1:length(diffsig2)-1)<0);
    if ~isempty(lastind)
        onsind = [onsind lastdur(lastind+1)];
    elseif length(ressig)-onsind(end)>150
        onsind = [onsind length(ressig)];
    end
end

freqdis = zeros(1,size(onsind)-1);

medwinsiz2 = 25;
ressig = detsig-medfilt1(detsig,medwinsiz2);

peaksig = ressig;
peaksig(find(ressig<peakthr))=0;
dicrsig = ressig;
dicrsig(find(abs(ressig)<dicrthr))=0;
peaind = zeros(1,size(onsind)-1);
dicrind = zeros(size(peaind));
secpind = zeros(size(peaind));
for i=1:length(onsind)-1
   sigindrange = onsind(i)+1:onsind(i+1)-1;
   difsig = diff(peaksig(sigindrange));
   canpind = find(difsig(1:length(difsig)-1).*difsig(2:length(difsig))<0 &difsig(1:length(difsig)-1)>0);
   subtrind = find(detsig(sigindrange(canpind+1))<mean(detsig(sigindrange)));
   if ~isempty(subtrind) & (sigindrange(canpind(subtrind)+1)-onsind(i)<80 | onsind(i+1)-sigindrange(canpind(subtrind)+1)<20)
       canpind = setdiff(canpind,canpind(subtrind));
   end
   if length(canpind)>1
       if length(canpind)~=2
           freqdis(i) = -2;
           peaind(i) = -1;
           dicrind(i)=-1;
           secpind(i)=-1;
       else
           peakcan = sigindrange(canpind(1)+1);
           secpcan = sigindrange(canpind(2)+1);
           dicrange = peakcan+1:secpcan-1;
           difdicrsig = diff(dicrsig(dicrange));
           dicrcanind = find(difdicrsig(1:length(difdicrsig)-1).*difdicrsig(2:length(difdicrsig))<0 & difdicrsig(1:length(difdicrsig)-1)<0);
           if isempty(dicrcanind)
               dicrcanind = find(difdicrsig(1:length(difdicrsig)-1).*difdicrsig(2:length(difdicrsig))==0 & difdicrsig(1:length(difdicrsig)-1)<0);
           end
           dicrcan = dicrange(dicrcanind(1)+1);
           dicrhei = detsig(dicrcan);
           secphei = detsig(secpcan);
           peakhei = detsig(peakcan);
           onsehei = detsig(onsind(i));
           if (secphei-dicrhei)> (peakhei-onsehei)/3 || dicrhei>peakhei
               freqdis(i) = -1;
               peaind(i) = -1;
               dicrind(i)=-1;
               secpind(i)=-1;
           else
               peaind(i) = peakcan;
               dicrind(i)=dicrcan;
               secpind(i)=secpcan;
               if dicrhei > 2*peakhei/3
                   freqdis(i) = 0;
               else
                   freqdis(i) = 1;
               end
           end
       end      
   else
       peaind(i) = sigindrange(canpind(1)+1);
       freqdis(i) = 1;
       dicrind(i)=-1;
       secpind(i)=-1;
   end
end
end