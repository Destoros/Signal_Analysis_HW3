close all
clear all
clc


id = 25;


%(A) Task 1 -----------------------------------------------------------
Task = 1;
fS = 1E6;
x = signalgenerator1(id, Task); %generates stationary signal


NDFT = 2^11;

%windows
h = ones(1,NDFT);
% h = hann(NDFT);
% h = kaiser(NDFT,0.5)

%normalize window
 h = h/(norm(h));

%plot window
% figure
%     plot(h)
%     title('window')

%pwelch(x,window,noverlap,nfft)
[pxx, f] = pwelch(x,h,NDFT/2,NDFT,fS);    %window.....Fenster
                                          %noverlap...Number Overlap; NDFT/2 = 50% overlap
                                          %nfft.......FFT mit n pins
                                   
%Bartlett method: segment the data and overlap them
%Welch method: segment the data with overlap = 50%
                                        
% cvv0 = sum(h.^2); NO need due to normalizing
% 
% %correct scaling
% pxx = pxx/cvv0; %divide by Signal Power


%plot PSD
figure
subplot(2,1,1)
    plot(f,pxx)
    title('Non-parametric PSD estimation')
    xlabel('f')
    ylabel('S_{xx}(e^{jw})')
    hold on
    
    %SEEN FROM PLOT:
    %Bandwitdhs from 0 - 200KHz and 300 - 500 KHz
    %at 1.748E5 Hz  and 2.251E5 Hz are two peaks (prolly from eigensignals)
        
%Total Signal Power

%via 2nd order moment
TSP_SecondOrderMoment = var(x) %var = E( (x - E(x))^2 )

%via PSD
delta_x = max(diff(f));
TSP_PSD = sum(pxx*delta_x)%Window h is normalized, no need to correct scaling
% Total_Signal_Power_4 = sum(pxx2)




p = 110; %model order
[a, e] = aryule(x, p);
b = 1;
f2 = linspace(0,fS/2,fS/2);

[FilterModel, w] = freqz(b,a,f2,fS); 
%Sxx = |H(e^jw)|^2 * sigma^2
pxx_PAR = abs(FilterModel).^2 * e; 

pxx_PAR_scaled = pxx_PAR/numel(FilterModel);

delta_x = max(diff(w));
TSP_PAR_PSD = sum(pxx_PAR_scaled*delta_x)



%figure
subplot(2,1,2)
    plot(f2,abs(pxx_PAR_scaled))
    title(['Parametric PSD estimation for p =' num2str(p)])
    xlabel('f')
    ylabel('S_{xx}(e^{jw})')
    hold off
    %SEEN FROM PLOT: 
    %Peak at eigensignal frequencies are lower than in the NON-PAR PSD
    %Probably overmodelled, hence the oscillation in the plot
    %y-scale is different than from the NON-PAR PSD
    
    %MODIFICATIONS: The scaling is incorrect and the model order is not
    %correct
    
    
    
% (B) Task 2 -----------------------------------------------------------
Task = 2;
x = signalgenerator1(id, Task); %generates sinusoidal signal
f_S = 1E6; 

%H(z) schreiben
r = 0.95; %1 < r << 0
w0 = pi/4; %0 < w0 < pi

%with AR processes we can model narrow peaks by placing a pole close to the unit circle (single frequency)!



z = tf('z');
Hz = z^2/((z-r*exp(1i*w0))*(z-r*exp(-1i*w0))); % AR process

[b, a] = tfdata(Hz, 'v'); % get coefficients from AR  process, 'v' outputs them as an array and not cell

figure
    freqz(b,a)
    title('Signal generated from the AR process')
    

    

%Generate own signal
A = 1;
f = 13.242E03;

t = 0:1/f_S:(10000-1)/f_S;

x1 = A*sin(2*pi*f*t);

%estimate coefficients of AR(2) 
p = 2;
a1 = aryule(x1,p);    

a1_angle = angle(pole(tf(1,a1))); %gives the angle in the unit circle
est_freq1_PAR_PSD = a1_angle*f_S/2/pi; %pi "entspricht" f_S/2 -> Calculate f
est_freq1_PAR_PSD
%negativ component due to negative frequency, can be elided

%estiamted freq. is close to true frequency


%Comparsion via NON-PAR PSD

[px1x1, w1] = pwelch(x1,h,[],[],f_S); %if noverlap and nfft is not set, pwelch uses the welch method with nfft = length(x1)
figure
    plot(w1,px1x1)
    title('Frequency via Non-parametric PSD ')
    ylabel('Sxx');
    xlabel('frequency')
    
    
[~, x_pmax] = max(px1x1);
est_freq1_NON_PAR_PSD = w1(x_pmax)

error_freq_PAR_PSD = est_freq1_PAR_PSD(est_freq1_PAR_PSD>=0) - f
est_freq1_NON_PAR_PSD = est_freq1_NON_PAR_PSD - f




%estmate freq of given signal
p = 2; %AR(2) known from the task sheet
a2 = aryule(x,p);

a2_angle = angle(pole(tf(1,a2))); %gives the angle in the unit circle
est_freq2_given_signal = a2_angle*f_S/2/pi %pi is equivalent to f_S/2 -> Calculate f

%Comparison via FFT
X = fftshift(fft(x));
f_fft = linspace(-f_S/2,f_S/2,length(X));
figure
    plot(f_fft, abs(X))
    ylabel('|X(e^{jw})|')
    xlabel('f')
    title('Check frequency via fft')
    
    %SEEN FROM PLOT: f = 5.646E4 Hz
    %in the same range as est_freq2, so the estimation via angle of poles
    %is probably valid implemented
    



%% (C) Task 3 ---------------------------------------------------

z = tf('z');
H = ((z - 0.55)*(z + 0.25)*(z - 1i*0.85)*(z  +1i*0.85))^-1;

p = 4; %model order
[b0,a0] = tfdata(H,'v'); 

% n_points_pool = [100, 500, 1000];
n_points_pool = 100;
n_simulations = 10000;

for n_points = n_points_pool %obsolete right now
    for ii = 1:n_simulations

        x = randn(n_points,1);
        y = filter(b0,a0,x);
        y = y(:);

        %Yule-Walker unbiased
        [rub, mub] = xcorr(y,'unbiased');

    %   form toeplitzmatrix
    %     Rub = [ rub[0], rub[-1], rub[-2], ... , rub[1-p]
    %             rub[1], rub[0],  rub[-1], ... ,
    %              ...

        rub = rub(mub >= 0 & mub <= p); %rhs needs values of rxx beginning from [-1] not [-0]
        Rub = toeplitz(conj(rub(1:p))); %conj bc Rub has the structure of conjugate toeplitz matrix
        aub = -Rub^-1 * rub(2:end); %- bc of convention during derivation
        aub=[1,aub']; %add leading 1 to coefficients of a





        %Yule-Walker unbiased
        [rb, mb] = xcorr(y,'biased');

        rb = rb(mb >= 0 & mb <= p); %rhs needs values of rxx beginning from [-1] not [-0]
        Rb = toeplitz(conj(rb(1:p))); %conj bc Rub has the structure of conjugate toeplitz matrix
        ab = -Rb^-1 * rb(2:end); %- bc of convention during derivation
        ab=[1,ab']; %add leading 1 to coefficients of a



        %Least Squares

        %create Matrix
        %Y1 = [y[0], y[-1], y[-2], ...
        %      y[1], y[0],  y[-1], ...
        %      ...

        N = length(y);
        y1=y(p+1:N); %first p entries are needed for the values of the Y-Matrix with negative index
        Y1=toeplitz(y(p:N-1),y(p:-1:1).'); %compact way to generate matrix in wanted form


        % compute the AR coffecients of Lsq
        alsq = -pinv(Y1) * y1;
        alsq = [1,alsq'];
    %     blsq = zeros(1,p);
    %     blsq(1) = 1; %NOT NEEDED

        %use inbuild function
        ar_inbuild = aryule(y,p); %for comparison


        error_ub(ii,:) = aub - a0;
        error_b(ii,:) = ab - a0;
        error_lsq(ii,:) = alsq - a0;    
        error_inbuild(ii,:) = ar_inbuild - a0;
    end

    mean_ub = mean(error_ub,1);
    var_ub = var(error_ub,[],1);%variance along the first dimension, which is the variance over different realizaitons
    mse_ub = mean(error_ub.^2,1);

    mean_b = mean(error_b,1);
    var_b = var(error_b,[],1);
    mse_b = mean(error_b.^2,1);

    mean_lsq = mean(error_lsq,1);
    var_lsq = var(error_lsq,[],1);
    mse_lsq = mean(error_lsq.^2,1);

    mean_inbuild = mean(error_inbuild,1); %ADDITIONAL
    var_inbuild = var(error_inbuild,[],1);
    mse_inbuild = mean(error_inbuild.^2,1);
    


    cols = {};
    for ii = 0:p
        cols = [cols, ['a' num2str(ii)]];
    end
    


 

    names = {'Error_Unbiased', ' Error_Biased', ' Error_Least_Square'};
    T_mean = table(mean_ub.', mean_b.', mean_lsq.','RowNames', cols, 'VariableNames', names);
    T_var = table(var_ub.', var_b.', var_lsq.','RowNames', cols, 'VariableNames', names);
    T_mse = table(mse_ub.', mse_b.', mse_lsq.','RowNames', cols, 'VariableNames', names);

    writetable(T_mean,'mean.xlsx') %to easily include in presentation
    writetable(T_var,'var.xlsx') %to easily include in presentation
    writetable(T_mse,'mse.xlsx') %to easily include in presentation
    
    
%CREATE UI TABLE DIRECTLY IN MATLAB
%%%%%%
    create_tables = 0; %ON/OFF
%%%%%

    if create_tables == 1    
        
        rows_mean = {'Error Unbiased, mean', ' Error Biased, mean', ' Error Least Square, mean', 'Error Inbuild, mean'};
        rows_var = {'Error Unbiased, variance', ' Error Biased, variance', ' Error Least Square, variance', 'Error Inbuild, variance'};
        rows_mse = {'Error Unbiased, mse', ' Error Biased, mse', ' Error Least Square, mse', 'Error Inbuild, mse'};

        table_mean = [mean_ub; mean_b; mean_lsq; mean_inbuild];
        table_var = [var_ub; var_b; var_lsq; var_inbuild];
        table_mse = [mse_ub; mse_b; mse_lsq; mse_inbuild];


        fig = uifigure('Name', 'mean', 'Position',[500 500 750 350]);
        uit = uitable(fig, 'Data', table_mean,'Position',[20 20 250 500]);
        uit.ColumnName = cols;
        uit.RowName = rows_mean;
        uit.Position = [20 20 710 310];

        fig = uifigure('Name', 'var', 'Position',[500 500 750 350]);
        uit = uitable(fig, 'Data', table_var,'Position',[20 20 250 500]);
        uit.ColumnName = cols;
        uit.RowName = rows_var;
        uit.Position = [20 20 710 310];

        fig = uifigure('Name', 'mse', 'Position',[500 500 750 350]);
        uit = uitable(fig, 'Data', table_mse, 'Position',[20 20 250 500]);
        uit.ColumnName = cols;
        uit.RowName = rows_mse;
        uit.Position = [20 20 710 310];
    end

end
%HOW CAN I CLOSE UIFIGURES?




