function [Sim_Times] = Create_Time_Matrix_Int(Chirps_Per_Frame,Tx_Start_Time,Frame_Time,Inter_Frame_Gap,EoSim,Chirp_Time)
i=1;
k = length(Inter_Frame_Gap);
Curr_Time = 0;
while(Curr_Time<EoSim)
    if(i==1)
        Sim_Times(i,1)=Tx_Start_Time;
    else
        Sim_Times(i,1) = Sim_Times(i-1,1) + Frame_Time + Inter_Frame_Gap(mod(i-1,k)==0*k+mod(i-1,k));
    end
    for j = 2:Chirps_Per_Frame
        Sim_Times(i,j) = Sim_Times(i,j-1) + Chirp_Time;
    end
    Curr_Time = Sim_Times(i,j);
    i = i+1;
end