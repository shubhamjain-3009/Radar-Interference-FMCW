clear all;
close all;
load("temp_file1.mat");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rf = 1;
Rc = 28;
range_res = (dsamp_rate*Ts*c)/(2*BW*N_sample);
Max_Range = (dsamp_rate*Ts*c)/(2*BW);
range_xaxis = (1:N_sample)*range_res;
u = -0.5:1/Angle_FFT_Len:0.5-1/Angle_FFT_Len; %% x-axis in world of sin(theta)*(sep/lambda)
angle_vals = asind((lambda/Srx)*u).';
wi = (kaiser(N_sample,19))';
wi = wi/sum(wi);        
Num_Chirp_Disp = 5;

tic
for frm = 1:N_f
    for chp = 1:Num_Chirp_Disp
        Return_Frame_Number = frm*1;
        Return_Chirp_Number = chp*(Chirps_Per_Frame/Num_Chirp_Disp);
        
        adcn_mat_ds = squeeze(adcn_mat_ds_totalsim(Return_Frame_Number, Return_Chirp_Number,:,:));
        RangeFFT_mat = fftshift(fft(adcn_mat_ds.*wi,N_sample,2));
        RangeFFTn = RangeFFT_mat(1,:);  %% Since this is only Range FFT, we pick out the first antenna waveform directly
        
        [pks,locns,~,p] = findpeaks(db(RangeFFTn(1,:)),range_xaxis);
        Trial = Peak_Finder(pks,locns,p,N_ref);
        [~,Trial_ind] = intersect(range_xaxis,Trial);

        %%%%%%%%%%%%%%%%% Angle Calculation DSP %%%%%%%%%%%%%%%%%%%%%%%%%%%
        AngleFFT_mat = zeros(N_ref,Angle_FFT_Len);
        for i = 1:length(Trial)
            AngleFFT_mat(i,:) = fftshift(fft(RangeFFT_mat(:,Trial_ind(i)),Angle_FFT_Len,1));
        end

        %%%%%%%%% Display Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(4)
        subplot(5,5,(frm-1)*5+chp); 
        %plot(range_xaxis,db(RangeFFTn(1,:)));
        plot(((-length(RangeFFTn(1,:))/2):length(RangeFFTn(1,:))/2-1),db(RangeFFTn(1,:)));
        grid on;
        title(sprintf('Frame %d Chirp %d',frm,Return_Chirp_Number));
        %axis([0 300 -200 -140]);
        figure(5)
        subplot(5,5,(frm-1)*5+chp);
        plot(angle_vals,abs(AngleFFT_mat));
        grid on;
    end
end
toc


%%
%%%%%%%%%% Detect Error Frame %%%%%%%%%%%%%%%%%%%

for frm = 1:2
    for chp = 1:Chirps_Per_Frame
        Return_Frame_Number = frm*1;
        Return_Chirp_Number = chp;
        adcn_mat_ds_filt = squeeze(adcn_mat_ds_totalsim(Return_Frame_Number, Return_Chirp_Number,:,:));
        RangeFFT_mat_filt = fftshift(fft(adcn_mat_ds_filt.*wi,N_sample,2));
        RangeFFTn_filt = RangeFFT_mat_filt(1,:);  %% Since this is only Range FFT, we pick out the first antenna waveform directly
        x = RangeFFTn_filt(0.5*length(RangeFFTn_filt)+1:end);
        d = RangeFFTn_filt(1:0.5*length(RangeFFTn));
        sumpow(frm,chp) = sum(abs(d).^2);
    end
end
figure(12)
plot(db(sumpow(1,:)))
hold on
%plot(db(sumpow(2,:)))
hold off
ylabel("Magnitude of power in  negative frequency bins (dB)")
xlabel("Chirp Number")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%
% %%%%%%%%%% LMS %%%%%%%%%%%%%%%%%%%
% 
% for frm = 1:2
%     for chp = 1:Chirps_Per_Frame
%         Return_Frame_Number = frm*1;
%         Return_Chirp_Number = chp;
%         
%         adcn_mat_ds_filt = squeeze(adcn_mat_ds_totalsim(Return_Frame_Number, Return_Chirp_Number,:,:));
%         RangeFFT_mat_filt = fftshift(fft(adcn_mat_ds_filt.*wi,N_sample,2));
%         RangeFFTn_filt = RangeFFT_mat_filt(1,:);  %% Since this is only Range FFT, we pick out the first antenna waveform directly
%         x = RangeFFTn_filt(0.5*length(RangeFFTn_filt)+1:end);
%         d = conj(fliplr(RangeFFTn_filt(1:0.5*length(RangeFFTn))));
%         if(chp==1)
%             x_tot = x;
%             d_tot = d;
%         else
%             x_tot = [x_tot x];
%             d_tot = [d_tot d];
%         end
%     end
% end
% 
% % figure(12)
% % plot(db(x_tot));
% % hold on
% % plot(db(d_tot),LineWidth=2);
% % hold off
% 
% l = 32;
% mu = 0.008;
%         
% lms = dsp.LMSFilter('Length',l,'StepSize',mu);
% [y,e] = lms(x_tot',d_tot');
% 
% y2 = x_tot'-d_tot';
% 
% figure(12)
% subplot(1,3,1)
% plot(db(x_tot));
% hold on
% plot(db(d_tot));
% hold off
% ylim([-220 -120])
% subplot(1,3,2)
% plot(db(y));
% %ylim([-220 -120])
% subplot(1,3,3)
% plot(db(y2));
% ylim([-220 -120])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%% Plot a Chirp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Return_Frame_Number = Rf;
Return_Chirp_Number = Rc;
adcn_mat_ds = squeeze(adcn_mat_ds_totalsim(Return_Frame_Number, Return_Chirp_Number,:,:));
RangeFFT_mat = fft(adcn_mat_ds.*wi,N_sample,2);
RangeFFTn = RangeFFT_mat(1,:);

figure(6)
subplot(1,3,1);
plot(range_xaxis,db(RangeFFTn(1,:)));
grid on;
subplot(1,3,2);
plot(real(adcn_mat_ds(1,:)));
subplot(1,3,3)
x2 = (adcn_mat_ds(1,:));
spectrogram(x2,[],[],[],dsamp_rate,'yaxis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inter_Frame_Gap = [Inter_Frame_Gap 0];
t_idle = t_start*1e06;
ramp_time = Ts*1e06;
k = K_ego*1e-12;
N_c = Chirps_Per_Frame;
IFG = Inter_Frame_Gap*1e06;
N_f = N_f;
t_total = (((t_idle+ramp_time)*N_c));
adc_start = adc_start*1e06;

x1 = 0:ramp_time;
x1 = k*x1;
x1 = [zeros(1,t_idle) x1];

x1_wvfm = x1;
x1_adc = [zeros(1,t_idle+adc_start) ones(1,length(x1_wvfm)-(t_idle+adc_start))];
x1_adc_wvfm = x1_adc;
for i=1:N_c
  x1_wvfm = [x1_wvfm x1];
  x1_adc_wvfm = [x1_adc_wvfm x1_adc];
end
flen = zeros(1,N_f+1);
x1_wvfm_frm_fin = [x1_wvfm zeros(1,IFG(1))];
%flen(2) = length(x1_wvfm_frm_fin);
for j = 1:N_f-1
    flen(j+1) = length(x1_wvfm_frm_fin);
    x1_wvfm_frm_fin = [x1_wvfm_frm_fin x1_wvfm zeros(1,IFG(j+1))];
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
    inter(i).x2 = 0:ramp_time_2(i);
    inter(i).x2 = k_2*inter(i).x2;
    inter(i).x2 = [zeros(1,t_idle_2(i)) inter(i).x2];
    inter(i).x2_wvfm = inter(i).x2;
    for j=1:N_c_2
        inter(i).x2_wvfm = [inter(i).x2_wvfm inter(i).x2];
    end
    inter(i).x2_wvfm_frm_fin = [inter(i).x2_wvfm zeros(1,IFG_2(i,1))];
%     for j = 1:N_f_2-1
%         inter(i).x2_wvfm_frm_fin = [inter(i).x2_wvfm_frm_fin inter(i).x2_wvfm zeros(1,IFG_2(i,j+1))];
%     end
    len = length(inter(i).x2_wvfm_frm_fin);
    q_1 = length(Inter_Frame_Gap_Int(i));
    q_2 = 1;
    while(len<flen(N_f))
        temp = IFG_2(mod(q_2-1,q_1)==0*q_1+mod(q_2-1,q_1))
        inter(i).x2_wvfm_frm_fin = [inter(i).x2_wvfm_frm_fin inter(i).x2_wvfm zeros(1,temp)];
        len = length(inter(i).x2_wvfm_frm_fin);
        q_2 = q_2+1;
    end
    inter(i).x2_wvfm_frm_fin = [zeros(1,tx_start_int(i)) inter(i).x2_wvfm_frm_fin];
end

figure(1)
plot(x1_wvfm_frm_fin,'b')
hold on
for i = 1:N_int
    plot(inter(i).x2_wvfm_frm_fin,'r')
end
hold off

Frame_Num = Return_Frame_Number;
Chirp_Num = Return_Chirp_Number;
Chirp_Start = flen(Frame_Num) + length(x1)*(Chirp_Num-1)+1;
Chirp_Length = length(x1);
ego_wvfm = x1_wvfm_frm_fin(Chirp_Start-0*Chirp_Length:Chirp_Start + 1*Chirp_Length);
for i = 1:N_int
    inter(i).int_wvfm = inter(i).x2_wvfm_frm_fin(Chirp_Start-0*Chirp_Length:Chirp_Start + 1*Chirp_Length);
end
figure(2)
subplot(2,2,[1 3])
plot((ego_wvfm),'b');
hold on
plot(ego_wvfm-15,'LineStyle',':','Color','black');
plot(ego_wvfm+15,'LineStyle',':','Color','black');
for i = 1:N_int
    plot(inter(i).int_wvfm,'r')
end
xline(adc_start+1);
xline(Chirp_Time*1e06+1);
hold off
subplot(2,2,2)
plot((ego_wvfm),'b');
hold on
plot(ego_wvfm-15,'LineStyle',':','Color','black');
plot(ego_wvfm+15,'LineStyle',':','Color','black');
for i = 1:N_int
    plot(inter(i).int_wvfm,'r')
end
axis([adc_start+1 Chirp_Time*1e06+1 0 350])
xline(adc_start+1);
xline(Chirp_Time*1e06+1);
hold off
title("Frame "+Rf+" ,Chirp "+Rc)
subplot(2,2,4)
plot(real(adcn_mat_ds(1,:)));
axis([0 512 -inf inf]);