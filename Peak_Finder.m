function [Locations] = Peak_Finder(pks,locns,p,numpeaks)
    
[p,idx] = sort(p,"descend");
locns = locns(idx);
pks = pks(idx);
for i=1:length(p)
    if(p(i)<1)
        break
    end
end
locns = locns(1:i);
pks = pks(1:i);
[~,idx] = sort(pks,"descend");
locns = locns(idx);
if (length(locns)>numpeaks)
    Locations = locns(1:numpeaks);
else
    Locations = locns(:);
end


