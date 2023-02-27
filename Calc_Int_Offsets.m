 function[s0,l0] = Calc_Int_Offsets(ToS,d,t_samp,Chirp_Time,Int_Time_Matrix,Tx_Start,Number_Of_Chirps)
global c;
ToS = ToS-(Tx_Start+d./c);
for i = 1:size(ToS,1)
    for j = 1:size(ToS,2)
        idx1 = find(Int_Time_Matrix(:,1)<ToS(i,j));
        if(isempty(idx1))
            lead_offset = abs(ToS(i,j));
            to=0;
        else
            idx2 = find(Int_Time_Matrix(idx1(end),:)<ToS(i,j));
            if(idx2(end)==Number_Of_Chirps && ToS(i,j)>(Int_Time_Matrix(idx1(end),idx2(end))+Chirp_Time))
                lead_offset = Int_Time_Matrix(idx1(end)+1,1)-ToS(i,j);
                to=0;
            else
                to = ToS(i,j)-Int_Time_Matrix(idx1(end),idx2(end));
                lead_offset = 0;
            end
        end
    end
end

l0 = floor(lead_offset/t_samp);
s0 = floor(to/t_samp);
