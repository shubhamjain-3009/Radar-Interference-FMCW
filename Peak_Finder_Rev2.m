function [Locations] = Peak_Finder_Rev2(pks,locns,p)
% [p,idx] = sort(p,"descend");
% locns = locns(idx);
% %pks = pks(idx);
% for i=1:length(p)
%     if(p(i)<1)  %%%%Threshold
%         break
%     end
% end
% Locations = locns(1:i);
%Peaks = pks(1:i);

% [~,idx] = sort(pks,"descend");
% locns = locns(idx);
% if (length(locns)>numpeaks)
%     Locations = locns(1:numpeaks);
% else
%     Locations = locns(:);
% end
k=1;
for i=1:length(p)
    if(p(i)>5)  %%%%Threshold
        Locations(k) = locns(i);
        k = k+1;
    end
end
