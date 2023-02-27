clear all;

load("Output_Para_1.mat");
range_res = (dsamp_rate*Ts*c)/(2*BW*N_sample);
Max_Range = (dsamp_rate*Ts*c)/(2*BW);
range_xaxis = (1:N_sample)*range_res;
u = -0.5:1/Angle_FFT_Len:0.5-1/Angle_FFT_Len; %% x-axis in world of sin(theta)*(sep/lambda)
angle_vals = asind((lambda/Srx)*u).';
        
Num_Chirp_Disp = 5;
tic
for frm = 1:N_f
    for chp = 1:Num_Chirp_Disp
        Return_Frame_Number = frm*1;
        Return_Chirp_Number = chp*(Chirps_Per_Frame/Num_Chirp_Disp);
        
        adcn_mat_ds = squeeze(adcn_mat_ds_totalsim(Return_Frame_Number, Return_Chirp_Number,:,:));
        RangeFFT_mat = fft(adcn_mat_ds,N_sample,2);
        RangeFFTn = RangeFFT_mat(1,:);
     
        [pks,locns,~,p] = findpeaks(db(RangeFFTn(1,:)),range_xaxis);
        Trial = Peak_Finder(pks,locns,p,N_ref);
        [~,Trial_ind] = intersect(range_xaxis,Trial);

        %%%%%%%%%%%%%%%%% Angle Calculation DSP %%%%%%%%%%%%%%%%%%%%%%%%%%%
        AngleFFT_mat = zeros(N_ref,Angle_FFT_Len);
        for i = 1:length(Trial)
            AngleFFT_mat(i,:) = fftshift(fft(RangeFFT_mat(:,Trial_ind(i)),Angle_FFT_Len,1));
        end

        %%%%%%%%% Display Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(2)
        subplot(5,5,(frm-1)*5+chp); 
        plot(range_xaxis,db(RangeFFTn(1,:)));
        grid on;
        %axis([0 300 -200 -140]);
        figure(3)
        subplot(5,5,(frm-1)*5+chp);
        plot(angle_vals,abs(AngleFFT_mat));
        grid on;
    end
end
toc



