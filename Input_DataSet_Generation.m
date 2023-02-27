function[Data_ID] = Input_DataSet_Generation(Num_Scenarios,Num_SubScenarios,Num_tar,Num_int,Override)
tic;
dmin = 5;
dmax = 220;
vmax = 35;
vmin = -35;
thetamax = 90;
thetamin = -90;
N_f = 5; %%% Number of frame for which the simulator needs to be run
Kmin = 8;
Kmax = 14;
f_start = 77e09;
Tx_Start_Time = 0;
K_ego = 10.71e12;
Noise_amp = 0.000000002;
margin_d = 5;
margin_theta = 5;
margin_v = 5;
margin_K = 0.001;

Data_ID = 1;

for i=1:Num_Scenarios
    if (Override == 0)
        N_tar = randi([3 10],1,1);
        N_int = randi([1 5],1,1);
    elseif(Override == 1)
        N_tar = Num_tar;
        N_int = Num_int;
    end
    for j=1:Num_SubScenarios
        
        
        rcs_int = ones(1,N_int);
        rcs_tar = ones(1,N_tar);
        Inter_Frame_Gap_Int_val = 20e-06*ones(1,N_int);
        Inter_Frame_Gap_Int = zeros(N_int,N_f-1);
        f_start_int = 77e09*ones(1,N_int);
        Tx_Start_Time_Int = 0*ones(1,N_int);


        d_int = floor(Rand_modification(dmin,dmax,N_int,margin_d));
        v_int = floor(Rand_modification(vmin,vmax,N_int,margin_v));
        theta_int = floor(Rand_modification(thetamin,thetamax,N_int,margin_theta));
        K_int = Rand_modification(Kmin,Kmax,N_int,margin_K)*1e12;        
%       d_int = dmin + (dmax-dmin).*rand(1,N_int); 
%       v_int = vmin + (vmax-vmin).*rand(1,N_int); 
%       theta_int = thetamin + (thetamax-thetamin).*rand(1,N_int); 
%       K_int = (Kmin + (Kmax-Kmin).*rand(N_int,1))*1e12;
        
        d_tar = floor(Rand_modification(dmin,dmax,N_int,margin_d));
        v_tar = floor(Rand_modification(vmin,vmax,N_int,margin_v));
        theta_tar = floor(Rand_modification(thetamin,thetamax,N_int,margin_theta));
%       d_tar = dmin + (dmax-dmin).*rand(1,N_tar); 
%       v_tar = vmin + (vmax-vmin).*rand(1,N_tar); 
%       theta_tar = thetamin + (thetamax-thetamin).*rand(1,N_tar);
        
        for l = 1:2
            if(l==1)
                int_present = 0;
                K_inter = K_ego*ones(1,N_int);
                Inter_Frame_Gap = 20e-06.*ones(1,N_f-1);
                save("Input_Para"+Data_ID+".mat");
                Data_ID = Data_ID+1;
            else
                int_present = 1;
                for m = 1:2
                    if(m==1) %%%same k
                        K_inter = K_ego*ones(1,N_int);
                        for q = 1:3
                            switch(q)
                               case 1
                                    Inter_Frame_Gap = 20e-06.*ones(1,N_f-1);
                                    save("Input_Para"+Data_ID+".mat");
                                    Data_ID = Data_ID+1;
                               case 2
                                    Inter_Frame_Gap = 1200e-06.*ones(1,N_f-1);
                                    save("Input_Para"+Data_ID+".mat");
                                    Data_ID = Data_ID+1;
                               case 3
                                    Inter_Frame_Gap = [800e-06 1200e-06 500e-06 600e-06];
                                    save("Input_Para"+Data_ID+".mat");
                                    Data_ID = Data_ID+1;
                            end
                        end
                    else %%%varying k
                        K_inter = K_int; 
                        for q = 1:3
                            switch(q)
                               case 1
                                    Inter_Frame_Gap = 20e-06.*ones(1,N_f-1);
                                     save("Input_Para"+Data_ID+".mat");
                                    Data_ID = Data_ID+1;
                               case 2
                                    Inter_Frame_Gap = 1200e-06.*ones(1,N_f-1);
                                     save("Input_Para"+Data_ID+".mat");
                                    Data_ID = Data_ID+1;
                               case 3
                                    Inter_Frame_Gap = [800e-06 1200e-06 500e-06 600e-06];
                                     save("Input_Para"+Data_ID+".mat");
                                    Data_ID = Data_ID+1;
                            end
                        end
                    end
                end
            end
        end 
    end
end
Time_Elapsed = toc;