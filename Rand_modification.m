% d_min = 2;
% d_max = 220;
% N_int = 5;
% 
% 
% min = d_min;
% max = d_max;
% n = N_int;
% d_int = y;

function [y] = Rand_modification(min,max,n,margin)
y = zeros(1,n);
for i = 1:n
    y(i) = min + (max-min).*rand(1,1);
    flag =1;
    while(flag==1 && i>1)
        flag=0;
        for j = 1:i-1
            if(abs(y(j)-y(i))<margin)
                flag=1;
            end
        end
        if (flag==1)
            y(i) = min + (max-min).*rand(1,1);
        end
    end
end

