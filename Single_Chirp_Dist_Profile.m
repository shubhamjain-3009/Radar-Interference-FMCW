clear all;
close all;
%%%%%%%%% Defining System Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global c;
c = 3e08;
Dec_factor = 1024; %%Digital to Analog Rate multiplication factor
dsamp_rate = 20e06; %%ADC Sampling Rate
asamp_rate = dsamp_rate*Dec_factor; %%Equivalent Analog Sampling Rate
t_samp = 1/asamp_rate; %%Sample time in the equivalent Analog Domain

%%%%%%%%% Defining Victim Radar Transmit Variables %%%%%%%%%%%%%%%%%%%%%%%%
f_start = 77e09;
Ts = 28e-06; %%ramp time. Total time = t_start+Ts
BW = 300e06;
k = BW/Ts;
%k = 10.71e12;
%BW = k*Ts;
t_start = 8e-06;
Chirps_Per_Frame=50;
Frame_Time = (Ts+t_start)*Chirps_Per_Frame;
Chirp_Time = (Ts+t_start);

%%%%%%%%% Defining Victim Radar Receive Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nrx = 16;
lambda = c/f_start;
Srx = lambda/2;
adc_start = 10e-06; %% must be greater than t_start
N_sample = 512; %% Ts@20M > 25.6us (N_sample*samp_time)
adc_len = N_sample/dsamp_rate;
Angle_FFT_Len = 1024;

%%%%%%%%% Defining Simulation Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = 1;
Time_of_Sim = Chirp_Time; %%ToS > (2*d/c)@ d=Max_Loc, => ToS>((samp_rate*Ts)/BW) 
att_factor_tar = 4; %%%% exponent of d for target. Typically d^4, but here it is d^2. If no att, then make it 0. 

%%%%%%%%%%%%%%%%%%%Transmit Waveform Creation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = (0:t_samp:Ts-t_samp);
x = exp(1j*2*pi*(f_start + 0.5*k.*t).*t);
x = [zeros(1,floor(t_start/t_samp)) x];
x_adc = x(floor(adc_start/t_samp):floor((adc_len+adc_start)/t_samp));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

range_res = (dsamp_rate*Ts*c)/(2*BW*N_sample);
range_xaxis = (1:N_sample)*range_res;
pk = zeros(1,N_sample);

for i = 1:N_sample
     
    %%%%%%%%%%%% Sample Offset Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    d = range_xaxis(i);
    ToS = Time_of_Sim-(2.*d/c);
    TStep = (Ts+t_start)*Chirps_Per_Frame;
    Coarse = TStep.*floor(ToS./TStep);
    ToS = ToS-Coarse;
    TStep = (Ts+t_start);
    Fine = TStep.*floor(ToS./TStep);
    to = ToS-Fine;
    
    %%%%%%%%%%%%%%%%% Receive Waveform Creation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s0 = floor(to./t_samp);
    Y_mat = zeros(1,length(x));
    
    xn = G.*([x((s0(1)+1):length(x)) x(1:(s0(1)))]/(d(1)^att_factor_tar));
    Y_mat(1,:) = xn;
    
    Y_Rx = Y_mat;

    %%%%%%%%%%%%% Downconversion (Mixer) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    Y_Rx = Y_Rx(floor(adc_start/t_samp):floor((adc_len+adc_start)/t_samp));
    adcn_mat = conj(Y_Rx).*(x_adc);
 
    %%%%%%%%%%%%%%%%%%% ADC Sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    adcn_mat_ds = decimate(adcn_mat,Dec_factor);
    adcn_mat_ds = adcn_mat_ds(:,1:N_sample);
    w = (kaiser(N_sample,19))';
    w = w/sum(w);
    adcn_mat_ds_w = adcn_mat_ds.*w;

    %%%%%%%%%%%%%%%%% Range Calculation DSP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RangeFFT_mat = fft(adcn_mat_ds_w,N_sample,2);
    RangeFFTn = RangeFFT_mat(1,:);
    [apk,aloc,~,~]=findpeaks(db(RangeFFTn(1,:)),range_xaxis,'SortStr','descend');
    loc = aloc(1);
    pk(i) = apk(1);
    i
    toc
end
Max_Range = (dsamp_rate*Ts*c)/(2*BW);
ulim = find(range_xaxis>220);
llim = find(range_xaxis>1);
plot(range_xaxis(llim(1):ulim(1)),pk(llim(1):ulim(1)));

