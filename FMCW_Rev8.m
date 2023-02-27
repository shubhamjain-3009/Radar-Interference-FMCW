clear all;
close all;
%%%%%%%%% Defining Simulation Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_f = 2; %%% Number of frame for which the simulator needs to be run

d_int = [75];
N_int = length(d_int);%%%%%Comment out hen turning into a function
v_int = [-22.1989];
rcs_int = ones(1,N_int);
theta_int = [52];
int_present = 1; %1 = Present, 0 = Not 
K_inter = [11e12];
f_start_int = 8e09*ones(1,N_int);
Inter_Frame_Gap_Int_val = 20e-06*ones(1,N_int);
Inter_Frame_Gap_Int = zeros(N_int,N_f-1);
Tx_Start_Time_Int = 0*ones(1,N_int);

d_tar = [110 40 28 220];
v_tar = [-5.3935 -7.9121 -17.3315 4];
theta_tar = [15 -30 -40 30];
rcs_tar = ones(1,length(d_tar));
K_ego = 11e12;
f_start = 8e09;
%Inter_Frame_Gap = 20e-06.*ones(1,N_f-1);
%Inter_Frame_Gap = 1200e-06.*ones(1,N_f-1);
Inter_Frame_Gap = [800e-06 1200e-06 500e-06 600e-06];
Tx_Start_Time = 0;

Noise_amp = 0.000000005;
%%%%%%%%%%%%%%%%%%% Sim Derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_int = length(d_int);
att_factor_Int = 2;
N_tar = length(d_tar);
att_factor_tar = 4; %%%% exponent of d for target
dist = [d_tar d_int]; %%Just add elements for multiple targets
v = [v_tar v_int];
theta_initial = [theta_tar theta_int];%% between +90 and -90
rcs = [rcs_tar rcs_int];
N_ref = length(dist);
tar_present = 1;

Noise_present = 1;

%%%%%%%%% Defining System Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
global c;
c = 3e08;
Dec_factor = 1024; %%Digital to Analog Rate multiplication factor
dsamp_rate = 20e06; %%ADC Sampling Rate
asamp_rate = dsamp_rate*Dec_factor; %%Equivalent Analog Sampling Rate
t_samp = 1/asamp_rate; %%Sample time in the equivalent Analog Domain

%%%%%%%%% Defining Victim Radar Receive Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nrx = 16; %%Number of Antennae
N_sample = 512; %% Ts@20M > 25.6us (N_sample*samp_time)
adc_len = N_sample/dsamp_rate; %% Length of ADC capture
Angle_FFT_Len = 1024;

%%%%%%%%% Defining Victim Radar Variables %%%%%%%%%%%%%%%%%%%%%%%%
Ts = 28e-06; %%ramp time. Total time = t_start+Ts 
BW = Ts*K_ego;
%Ts = BW/K_ego;
lambda = c/f_start;
Srx = lambda/2;
adc_start = 10e-06; %% must be greater than t_start
%k = BW/Ts;
k = K_ego;
%k = 10.71e12;
%BW = k*Ts;
t_start = 8e-06;
Chirps_Per_Frame=50;
Chirp_Time = (Ts+t_start);
Frame_Time = Chirp_Time*Chirps_Per_Frame;
Sim_Times = Create_Time_Matrix(N_f,Chirps_Per_Frame,Tx_Start_Time,Frame_Time,Inter_Frame_Gap,Chirp_Time);
EoSim = Sim_Times(end,end)+Chirp_Time;

%%%%%%%%%%%%%%%%%%%%%%% Scenario Generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = zeros(N_f,Chirps_Per_Frame,N_ref);
theta = zeros(N_f,Chirps_Per_Frame,N_ref);
vf = zeros(N_f,Chirps_Per_Frame,N_ref);
for i=1:size(Sim_Times,1)
    for j=1:size(Sim_Times,2)
        ToS = Sim_Times(i,j);
        [d(i,j,:),theta(i,j,:),vf(i,j,:)] = Update_Scene_mod(dist,theta_initial,v,ToS);
    end
end

A = zeros(N_f,Chirps_Per_Frame,Nrx,N_ref);
mu = (-2*pi*f_start*Srx.*sind(theta))/c ;%% Needs to be made robust for different f_starts of interferer and victim
for i = 1:Nrx
    A(:,:,i,:) = exp(1j*mu*(i-1)); 
end
A (:,:,:,(N_ref+1:N_ref+N_int)) = A (:,:,:,(N_ref-N_int+1:N_ref));

%%%%%%%%% Defining Interferer Radar Transmit Variables %%%%%%%%%%%%%%%%%%%%%
Ts_int = 28e-06*ones(1,N_int); %%% TS 20 to 30
BW_int = Ts_int.*K_inter;
BW_int = 300e06*ones(1,N_int);
Ts_int = (BW_int./K_inter); %%% TS 20 to 30
%k_int = BW_int./Ts_int;
k_int = K_inter;
t_start_int = 8e-06*ones(1,N_int);
Chirps_Per_Frame_int = 50*ones(1,N_int);
for i=1:N_int
    Inter_Frame_Gap_Int(i,:) = Inter_Frame_Gap_Int_val(i).*ones(1,N_f-1);
end
Chirp_Time_int = (Ts_int+t_start_int);
Frame_Time_int = Chirp_Time_int.*Chirps_Per_Frame_int;

%%%%%%%%%%%%%%% Interferer Transmit Waveform Creation %%%%%%%%%%%%%%%%%%%%%
Int_Offset_s = zeros(N_f,Chirps_Per_Frame,N_int);
Int_Offset_l = zeros(N_f,Chirps_Per_Frame,N_int);
inter = struct();
for i=1:N_int
    t = (0:t_samp:Ts_int(i)-t_samp);
    x_int = exp(1j*2*pi*(f_start_int(i) + 0.5*k_int(i).*t).*t);
    x_int = [zeros(1,floor(t_start_int(i)/t_samp)) x_int];
    inter(i).y_int = [x_int x_int x_int];
    inter(i).Int_Time_Matrix(:,:) = Create_Time_Matrix_Int(Chirps_Per_Frame_int(i),Tx_Start_Time_Int(i),Frame_Time_int(i),Inter_Frame_Gap_Int(i,:),EoSim,Chirp_Time_int(i));
end

%%%%%%%%%%%%%%%%%% Ego Vehicle Transmit Waveform Creation %%%%%%%%%%%%%%%%%
t = (0:t_samp:Ts-t_samp);
%x = exp(1j*2*pi*(f_start + 0.5*k.*t).*t);
x = cos(2*pi*(f_start + 0.5*k.*t).*t);
x = [zeros(1,floor(t_start/t_samp)) x];
len = length(x);
%x_adc = x(floor(adc_start/t_samp):floor((adc_len+adc_start)/t_samp));

xr = exp(1j*2*pi*(f_start + 0.5*k.*t).*t);
xr = [zeros(1,floor(t_start/t_samp)) xr];
%xr = x;
x_adc = xr(floor(adc_start/t_samp):floor((adc_len+adc_start)/t_samp));

%%%%%%%%%%%%%%%%%% Time Offset Array Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%
adcn_mat_ds_frame = zeros(Chirps_Per_Frame,Nrx,N_sample);
%adcn_mat_SNR_ds_frame = zeros(Chirps_Per_Frame,Nrx,N_sample);
adcn_mat_noint_ds_frame = zeros(Chirps_Per_Frame,Nrx,N_sample);

adcn_mat_ds_totalsim = zeros(N_f,Chirps_Per_Frame,Nrx,N_sample);
%adcn_mat_SNR_ds_totalsim = zeros(N_f,Chirps_Per_Frame,Nrx,N_sample);
adcn_mat_noint_ds_totalsim = zeros(N_f,Chirps_Per_Frame,Nrx,N_sample);

Ref_Offset_s = zeros(N_f,Chirps_Per_Frame,N_ref);
Ref_Offset_l = zeros(N_f,Chirps_Per_Frame,N_ref);
lead_offset = zeros(1,N_ref);
to = zeros(1,N_ref);

for i = 1:N_f
    for j=1:Chirps_Per_Frame
        for l = 1:N_int
        [Int_Offset_s(i,j,l),Int_Offset_l(i,j,l)] = Calc_Int_Offsets(Sim_Times(i,j),d(i,j,end-N_int+l),t_samp,Chirp_Time_int(l),inter(l).Int_Time_Matrix,Tx_Start_Time_Int(l),Chirps_Per_Frame_int(l));
        end
        for l = 1:N_ref
            if(Sim_Times(i,j)<2*d(i,j,l)/c)
                lead_offset(l) = (2*d(i,j,l)/c)-Sim_Times(i,j);
                to(l)=0;
            else
                lead_offset(l) = 0;
                ToS = Sim_Times(i,j)-(2*d(i,j,l)/c)-Sim_Times(i,1);
                to(l) = ToS-Chirp_Time*floor(ToS/Chirp_Time);
            end
        end
        Ref_Offset_s(i,j,:) = floor(to./t_samp);
        Ref_Offset_l(i,j,:) = floor(lead_offset./t_samp);
    end
end
toc

%%%%%%%%%%%%%%%%%%% Main Receiver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_exec = zeros(1,N_f);
Interference = zeros(N_int,len);
for j = 1:N_f
    tic
    for m = 1:Chirps_Per_Frame
        %%%%%%%%%%%% Sample Offset Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s0 = Ref_Offset_s(j,m,:);
        l0 = Ref_Offset_l(j,m,:);
        %%%%%%%%%%%%%%%%% Receive Waveform Creation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Y_mat = zeros(N_ref+N_int,len);
        %Y_mat_SNR = zeros(N_ref+N_int,len);

        for i=1:N_ref
            if(l0(i)==0)
                xn = vf(j,m,i).*([x((s0(i)+1):len) x(1:(s0(i)))]/(d(j,m,i)^att_factor_tar));
            else
                if(l0(i)>=len)
                    xn = zeros(1,len);
                else
                    xn = [zeros(1,l0(i)) x(1:(len-l0(i)))/(d(j,m,i)^att_factor_tar)];
                end
            end
            Y_mat(i,:) = rcs(i)*tar_present*xn;%%%  control for targets
        end
        
        for i = 1:N_int
            if(Int_Offset_l(j,m,i)==0)
                xn_int = inter(i).y_int(1,Int_Offset_s(j,m,i)+1:Int_Offset_s(j,m,i)+len);
            else
                if(Int_Offset_l(j,m,i)>=len)
                    xn_int = zeros(1,len);
                else
                    xn_int = [zeros(1,Int_Offset_l(j,m,i)) inter(i).y_int(1,1:(len-Int_Offset_l(j,m,i)))];
                end
            end
            Interference(i,:) = vf(j,m,end)*xn_int/(d(j,m,end)^att_factor_Int);
        end

        Y_mat_noint = Y_mat;
        Y_mat(N_ref+1:end,:) = int_present*Interference(:,:); %%%  cotrol for interferer
        %Y_mat_SNR(N_ref+1:end,:) = int_present*Interference(:,:); %%%  cotrol for interferer
        
        Noise_wvfm = Noise_amp*(2*rand(N_ref+N_int,len)-1+2j*rand(N_ref+N_int,len)-1j);
        
        Y_mat = Y_mat + Noise_present*Noise_wvfm;%%%%%%% control for noise
        %Y_mat_SNR = Y_mat_SNR + Noise_present*Noise_wvfm;%%%%%%% control for noise
        Y_mat_noint = Y_mat_noint + Noise_present*Noise_wvfm;
        
        Y_Rx = squeeze(A(j,m,:,:))*Y_mat;
        %Y_Rx_SNR = squeeze(A(j,m,:,:))*Y_mat_SNR;
        Y_Rx_noint = squeeze(A(j,m,:,:))*Y_mat_noint;
       
        %%%%%%%%%%%%% Downconversion (Mixer) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Y_Rx = Y_Rx(:,floor(adc_start/t_samp):floor((adc_len+adc_start)/t_samp));
        %Y_Rx_SNR = Y_Rx_SNR(:,floor(adc_start/t_samp):floor((adc_len+adc_start)/t_samp));
        Y_Rx_noint = Y_Rx_noint(:,floor(adc_start/t_samp):floor((adc_len+adc_start)/t_samp));

        adcn_mat = conj(Y_Rx).*(x_adc);
        %adcn_mat_SNR = conj(Y_Rx_SNR).*(x_adc);
        adcn_mat_noint = conj(Y_Rx_noint).*(x_adc);
        %%%%%%%%%%%%%%%%%%% ADC Sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        adcn_mat_ds = zeros(Nrx,N_sample+1);
        %adcn_mat_SNR_ds = zeros(Nrx,N_sample+1);
        adcn_mat_noint_ds = zeros(Nrx,N_sample+1);
        for i=1:Nrx
        adcn_mat_ds(i,:) = decimate(adcn_mat(i,:),Dec_factor);
        %adcn_mat_SNR_ds(i,:) = decimate(adcn_mat_SNR(i,:),Dec_factor);
        adcn_mat_noint_ds(i,:) = decimate(adcn_mat_noint(i,:),Dec_factor);
        %adcn_mat_ds(i,:) = decimate(adcn_mat(i,:),Dec_factor,'fir');
        end
    adcn_mat_ds_frame(m,:,:)= adcn_mat_ds(:,1:N_sample);
    %adcn_mat_SNR_ds_frame(m,:,:)= adcn_mat_SNR_ds(:,1:N_sample);
    adcn_mat_noint_ds_frame(m,:,:) = adcn_mat_noint_ds(:,1:N_sample);
    end
adcn_mat_ds_totalsim(j,:,:,:) = adcn_mat_ds_frame;
%adcn_mat_SNR_ds_totalsim(j,:,:,:) = adcn_mat_SNR_ds_frame;
adcn_mat_noint_ds_totalsim(j,:,:,:) = adcn_mat_noint_ds_frame;

time_exec(j) = toc
end
%save("sim_wvfm_&_para.mat",'dsamp_rate','Ts','c',"BW","N_sample","Angle_FFT_Len","lambda","Srx","N_f","Chirps_Per_Frame","adcn_mat_SNR_ds_totalsim" ...
%    ,"N_ref","N_int","d","theta","Chirp_Time","adcn_mat_ds_totalsim","Sim_Times");
Avg_Sim_Time_Per_Frame = sum(time_exec)/N_f
save("temp_file1.mat");

 