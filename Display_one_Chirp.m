clear all;
%close all;
load("Output_Para_1.mat");
%%%%%%%%%%%%%% Plot a Chirp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
range_res = (dsamp_rate*Ts*c)/(2*BW*N_sample);
Max_Range = (dsamp_rate*Ts*c)/(2*BW);
range_xaxis = (1:N_sample)*range_res;
u = -0.5:1/Angle_FFT_Len:0.5-1/Angle_FFT_Len; %% x-axis in world of sin(theta)*(sep/lambda)
angle_vals = asind((lambda/Srx)*u).';
adc_len = N_sample/dsamp_rate;

Return_Frame_Number = 4;
Return_Chirp_Number = 42;
adcn_mat_ds = squeeze(adcn_mat_ds_totalsim(Return_Frame_Number, Return_Chirp_Number,:,:));
RangeFFT_mat = fft(adcn_mat_ds,N_sample,2);
RangeFFTn = RangeFFT_mat(1,:);

figure(4)
subplot(1,2,1);
plot(range_xaxis,db(RangeFFTn(1,:)));
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter_Frame_Gap_Int_val = 20e-06;
% x_inv = 0:adc_len*10^6;
% k_inv = k*10^-12;
% k_inv_int = k_int*10^-12;
% x_inv_int = 0:Ts_int*10^6;
% 
% time = Sim_Times(Return_Frame_Number,Return_Chirp_Number);
% time_t = floor(time/(Frame_Time_int+Inter_Frame_Gap_Int_val))*(Frame_Time_int+Inter_Frame_Gap_Int_val);
% time = time-time_t;
% time_t = floor(time/Chirp_Time_int)*(Chirp_Time_int);
% time = ((time - time_t)+adc_start)*10^6;
% 
% y1 = [zeros(1,t_start_int*10^6-1) k_inv_int*(x_inv_int)];
% y1 = [y1 y1 y1];
% y1 = y1(time:time+adc_len*10^6);
% y2 = k_inv*x_inv+(adc_start-t_start)*10^6*k_inv;
% 
% subplot(1,2,2)
% plot(x_inv,y2);
% grid on
% hold on
% plot(x_inv,y1,"o");
% plot(x_inv,y2+16,'black');
% plot(x_inv,y2-16,'black');
% hold off

figure(5)
findpeaks(db(RangeFFTn(1,:)),range_xaxis,'Annotate','extents')