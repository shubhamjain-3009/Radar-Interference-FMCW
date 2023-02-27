clear all;
close all;
load("temp_file1.mat");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rf = 2;
Rc = 9;

range_res = (dsamp_rate*Ts*c)/(2*BW*N_sample);
Max_Range = (dsamp_rate*Ts*c)/(2*BW);
range_xaxis = (1:N_sample)*range_res;
wi = (kaiser(N_sample,19))';
wi = wi/sum(wi); 

Return_Frame_Number = Rf;
Return_Chirp_Number = Rc;
adcn_mat_ds = squeeze(adcn_mat_ds_totalsim(Return_Frame_Number, Return_Chirp_Number,:,:));
adcn_mat_noint_ds = squeeze(adcn_mat_noint_ds_totalsim(Return_Frame_Number, Return_Chirp_Number,:,:));

RangeFFT_mat = fft(adcn_mat_ds.*wi,N_sample,2);
RangeFFT_mat_noint = fft(adcn_mat_noint_ds.*wi,N_sample,2);

RangeFFTn = RangeFFT_mat(1,:);
RangeFFTn_noint = RangeFFT_mat_noint(1,:);

figure(1)
subplot(2,2,1);
plot(range_xaxis,db(RangeFFTn(1,:)));
grid on;
subplot(2,2,2);
plot(real(adcn_mat_ds(1,:)));
subplot(2,2,3);
plot(range_xaxis,db(RangeFFTn_noint(1,:)));
grid on;
subplot(2,2,4);
plot(real(adcn_mat_noint_ds(1,:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inter_Frame_Gap = [Inter_Frame_Gap 0];
t_idle = t_start*1e06;
ramp_time = Ts*1e06;
k = K_ego*1e-12;
N_c = Chirps_Per_Frame;
IFG = [Inter_Frame_Gap 0]*1e06;
N_f = N_f;
t_total = (((t_idle+ramp_time)*N_c));
adc_start_disp = adc_start*1e06;

x1 = 1:ramp_time;
x1 = k*x1;
x1 = [zeros(1,t_idle) x1];

x1_wvfm = x1;
%x1_adc = [zeros(1,t_idle+adc_start_disp) ones(1,length(x1_wvfm)-(t_idle+adc_start_disp))];
%x1_adc_wvfm = x1_adc;
for i=2:N_c
  x1_wvfm = [x1_wvfm x1];
  %x1_adc_wvfm = [x1_adc_wvfm x1_adc];
end
flen = zeros(1,N_f+1);
x1_wvfm_frm_fin = [x1_wvfm zeros(1,IFG(1))];
flen(2) = length(x1_wvfm_frm_fin);
for j = 1:N_f-1
    x1_wvfm_frm_fin = [x1_wvfm_frm_fin x1_wvfm zeros(1,IFG(j+1))];
    flen(j+2) = length(x1_wvfm_frm_fin);
end

t_idle_2 = t_start_int*1e06;
ramp_time_2 = Ts_int*1e06;
k_2 = K_inter*1e-12;
N_c_2 = Chirps_Per_Frame_int;
IFG_2 = [Inter_Frame_Gap_Int zeros(N_int,1)]*1e06;
%IFG_2 = Inter_Frame_Gap_Int*1e06;
N_f_2 = N_f;
tx_start_int = Tx_Start_Time_Int*1e06; %%%relative to a frame

for i = 1:N_int
    inter(i).x2 = 1:ramp_time_2(i);
    inter(i).x2 = k_2*inter(i).x2;
    inter(i).x2 = [zeros(1,t_idle_2(i)) inter(i).x2];
    inter(i).x2_wvfm = inter(i).x2;
    for j=2:N_c_2
        inter(i).x2_wvfm = [inter(i).x2_wvfm inter(i).x2];
    end
    inter(i).x2_wvfm_frm_fin = [inter(i).x2_wvfm zeros(1,IFG_2(i,1))];
%     for j = 1:N_f_2-1
%         inter(i).x2_wvfm_frm_fin = [inter(i).x2_wvfm_frm_fin inter(i).x2_wvfm zeros(1,IFG_2(i,j+1))];
%     end
    len = length(inter(i).x2_wvfm_frm_fin);
    q_1 = length(Inter_Frame_Gap_Int(i));
    q_2 = 1;
    while(len<flen(N_f+1))
        temp = IFG_2(mod(q_2-1,q_1)==0*q_1+mod(q_2-1,q_1));
        inter(i).x2_wvfm_frm_fin = [inter(i).x2_wvfm_frm_fin inter(i).x2_wvfm zeros(1,temp)];
        len = length(inter(i).x2_wvfm_frm_fin);
        q_2 = q_2+1;
    end
    inter(i).x2_wvfm_frm_fin = [zeros(1,tx_start_int(i)) inter(i).x2_wvfm_frm_fin];
end

figure(2)
plot(x1_wvfm_frm_fin,'b')
hold on
for i = 1:N_int
    plot(inter(i).x2_wvfm_frm_fin,'r')
end
hold off

Frame_Num = Return_Frame_Number;
Chirp_Num = Return_Chirp_Number;
Chirp_Start = flen(Frame_Num) + length(x1)*(Chirp_Num-1)+1
Chirp_Length = length(x1)
ego_wvfm = x1_wvfm_frm_fin(Chirp_Start-0*Chirp_Length:Chirp_Start + 1*Chirp_Length);
for i = 1:N_int
    inter(i).int_wvfm = inter(i).x2_wvfm_frm_fin(Chirp_Start-0*Chirp_Length:Chirp_Start + 1*Chirp_Length);
end
figure(3)
subplot(3,2,[1 3 5])
plot((ego_wvfm),'b');
hold on
plot(ego_wvfm-15,'LineStyle',':','Color','black');
plot(ego_wvfm+15,'LineStyle',':','Color','black');
for i = 1:N_int
    plot(inter(i).int_wvfm,'r')
end
xline(adc_start_disp+1);
xline(Chirp_Time*1e06+1);
hold off
subplot(3,2,2)
plot((ego_wvfm),'b');
hold on
plot(ego_wvfm-15,'LineStyle',':','Color','black');
plot(ego_wvfm+15,'LineStyle',':','Color','black');
for i = 1:N_int
    plot(inter(i).int_wvfm,'r')
end
axis([adc_start_disp+1 Chirp_Time*1e06 -inf inf])
hold off
title("Frame "+Rf+" ,Chirp "+Rc)
subplot(3,2,4)
plot(real(adcn_mat_ds(1,:)));
axis([1 512 -inf inf]);
subplot(3,2,6)
plot(real(adcn_mat_noint_ds(1,:)));
axis([0 512 -inf inf]);

