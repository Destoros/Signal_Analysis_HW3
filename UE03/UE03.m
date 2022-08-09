% ue3

clear all, close all, clc
disp('########### init ###########')

load Data_Ex3.mat
s2 = s2';
s2n = s2n';

% Set wavelet name. 
wname = {'haar','db10','sym8','bior2.2','rbio2.2'};

%% b) Analyzing and reconstructing
disp('########### b) Analyzing and reconstructing ###########')
% Compute the four filters associated with wavelet name given 
% by the input character vector wname. 
[Lo_D_h,Hi_D_h,Lo_R_h,Hi_R_h] = wfilters(wname{1});
[Lo_D_db,Hi_D_db,Lo_R_db,Hi_R_db] = wfilters(wname{2});
[Lo_D_sym,Hi_D_sym,Lo_R_sym,Hi_R_sym] = wfilters(wname{3});
[Lo_D_bi,Hi_D_bi,Lo_R_bi,Hi_R_bi] = wfilters(wname{4});
[Lo_D_rb,Hi_D_rb,Lo_R_rb,Hi_R_rb] = wfilters(wname{5});

% Haar-Wavelet
disp('########### Haar-filter ###########')
% deconstruction & reconstruction for lowpass filter path
h_low_D = upfirdn(s2,Lo_D_h,1,2); % no upsample before, downsample afterwards
h_low_R = upfirdn(h_low_D,Lo_R_h,2,1); % upsample before, no downsample afterwards

% deconstruction & reconstruction for highpass filter path
h_high_D = upfirdn(s1,Hi_D_h,1,2);
h_high_R = upfirdn(h_high_D,Hi_R_h,2,1);

test_x = h_low_R + h_high_R; % add lowpass and highpass path up

% delay between reconstructed and original
del = finddelay(s2, test_x);
string = sprintf('delay = %d', del);
disp(string)

figure, hold on, set(gca,'FontSize',26),set(gcf,'Color','White');
plot(test_x, 'LineWidth', 1.4)
grid on
plot(s2, 'LineWidth', 1.4);
legend('Reconstructed','Original')%, set(gca,'FontSize',16);
title('Reconstructing filter bank: Haar filter')%, set(gca,'FontSize',22);
axis tight

% db10-Wavelet
disp(['###########' wname{2} '-filter ###########'])
% deconstruction & reconstruction for lowpass filter path
h_low_D = upfirdn(s2,Lo_D_db,1,2); % no upsample before, downsample afterwards
h_low_R = upfirdn(h_low_D,Lo_R_db,2,1); % upsample before, no downsample afterwards

% deconstruction & reconstruction for highpass filter path
h_high_D = upfirdn(s1,Hi_D_db,1,2);
h_high_R = upfirdn(h_high_D,Hi_R_db,2,1);

test_x = h_low_R + h_high_R;

% delay between reconstructed and original
del = finddelay(s2, test_x);
string = sprintf('delay = %d', del);
disp(string)

figure, hold on, set(gca,'FontSize',26),set(gcf,'Color','White');
plot(test_x, 'LineWidth', 1.4)
grid on
plot(s2, 'LineWidth', 1.4);
legend('Reconstructed','Original')%, set(gca,'FontSize',16);
title(['Reconstructing filter bank:' wname{2} '  filter'])%, set(gca,'FontSize',22);
axis tight

% sym8-Wavelet
disp('########### sym8-filter ###########')
% deconstruction & reconstruction for lowpass filter path
h_low_D = upfirdn(s2,Lo_D_sym,1,2); % no upsample before, downsample afterwards
h_low_R = upfirdn(h_low_D,Lo_R_sym,2,1); % upsample before, no downsample afterwards

% deconstruction & reconstruction for highpass filter path
h_high_D = upfirdn(s1,Hi_D_sym,1,2);
h_high_R = upfirdn(h_high_D,Hi_R_sym,2,1);

test_x = h_low_R + h_high_R;

% delay between reconstructed and original
del = finddelay(s2, test_x);
string = sprintf('delay = %d', del);
disp(string)

figure, hold on, set(gca,'FontSize',26),set(gcf,'Color','White');
plot(test_x, 'LineWidth', 1.4)
grid on
plot(s2, 'LineWidth', 1.4);
legend('Reconstructed','Original')%, set(gca,'FontSize',16);
title('Reconstructing filter bank: sym8 filter')%, set(gca,'FontSize',22);
axis tight

% bior2.2-Wavelet
disp('########### bior2.2-filter ###########')
% deconstruction & reconstruction for lowpass filter path
h_low_D = upfirdn(s2,Lo_D_bi,1,2); % no upsample before, downsample afterwards
h_low_R = upfirdn(h_low_D,Lo_R_bi,2,1); % upsample before, no downsample afterwards

% deconstruction & reconstruction for highpass filter path
h_high_D = upfirdn(s1,Hi_D_bi,1,2);
h_high_R = upfirdn(h_high_D,Hi_R_bi,2,1);

test_x = h_low_R + h_high_R;

% delay between reconstructed and original
del = finddelay(s2, test_x);
string = sprintf('delay = %d', del);
disp(string)

figure, hold on, set(gca,'FontSize',26),set(gcf,'Color','White');
plot(test_x, 'LineWidth', 1.4)
grid on
plot(s2, 'LineWidth', 1.4);
legend('Reconstructed','Original')%, set(gca,'FontSize',16);
title('Reconstructing filter bank: bior2.2 filter')%, set(gca,'FontSize',22);
axis tight

% rbio2.2-Wavelet
disp('########### rbio2.2-filter ###########')
% deconstruction & reconstruction for lowpass filter path
h_low_D = upfirdn(s2,Lo_D_rb,1,2); % no upsample before, downsample afterwards
h_low_R = upfirdn(h_low_D,Lo_R_rb,2,1); % upsample before, no downsample afterwards

% deconstruction & reconstruction for highpass filter path
h_high_D = upfirdn(s1,Hi_D_rb,1,2);
h_high_R = upfirdn(h_high_D,Hi_R_rb,2,1);

test_x = h_low_R + h_high_R;

% delay between reconstructed and original
del = finddelay(s2, test_x);
string = sprintf('delay = %d', del);
disp(string)

figure, hold on, set(gca,'FontSize',26),set(gcf,'Color','White');
plot(test_x, 'LineWidth', 1.4)
hold on, grid on
plot(s2, 'LineWidth', 1.4);
legend('Reconstructed','Original')%, set(gca,'FontSize',16);
title('Reconstructing filter bank: rbio2.2 filter')%, set(gca,'FontSize',22);
axis tight


%% c) Denoising

disp('########### c) Denoising ###########')
wname_test = {'db45','db10','sym12', 'sym20'};
% wname = 'db8';
for jj = 1:length(wname_test)
    wname = wname_test{jj};
    [C,L] = wavedec(s2,5,wname);
    [Cn,Ln] = wavedec(s2n,5, wname);
    
    [rss, mss] = xcorr(s2);%S
    [rssn, mssn] = xcorr(s2n); %S+N    
    
    N = rssn(mssn==0) - rss(mss==0); %S+N - N;
    SNR_Input = 10*log10(rss(mss==0)/N)
    a = snr(s2n)
    %[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters('haar');

%     figure
%     plot(C)
%     hold on
%     plot(Cn)

    % highest highpass coefficient of noisy
    HP_L1 = Cn(Ln(end-1)+1:end);
    % second highest highpass coefficient of noisy
    HP_L2 = Cn(Ln(end-2)+1:Ln(end-1));

    hi_coeff = [HP_L2, HP_L1];

    eps = 3;

    hi_coeff = Cn(Ln(1)+Ln(2)+1:end);
    % soft tresholding
    Ln_cum = cumsum(Ln+1);
    Ln_cum(end) = [];
    start = Ln_cum(3);
    hard_threshold = 1;
    
    
    if hard_threshold == 1
        gamma = 0;        
    else
        gamma = eps;      
    end
    
    for i=start:numel(Cn)
        if(Cn(i) > eps)
            Cn(i) = Cn(i)-gamma;
        elseif(Cn(i) < -eps)
            Cn(i) = Cn(i)+gamma;
        else
            Cn(i) = 0;
        end
    end
    

    % figure
    % plot(hi_coeff)
    % hold on
    % plot(hi_coeff_new)
    % 
    % Cn_new = [Cn(1:Ln(end-2)), hi_coeff_new];
    % Cn_new = [Cn(1:(Ln(1)+Ln(2))), hi_coeff_new];
    S1n_rec = waverec(Cn, Ln, wname);
    

    [rSS, mSS] = xcorr(S1n_rec); %S+N    
    %S+N = rSS(mSS==0)
    %S = rSS(mSS==1) ?????????????????????????????????????????????????????
    
%     figure
%         plot(mSS,rSS)
        
    N2 = rSS(mSS==0) - rSS(mSS==1); %S+N - N;
    SNR_Output = 10*log10(rSS(mSS==1)/N2)
    b = snr(S1n_rec)
    

    figure, hold on, set(gca,'FontSize',26),set(gcf,'Color','White');
    plot(s2, 'LineWidth', 1.4)
    title(['Denoised Signal with ' wname ' filter'])
    plot(S1n_rec, 'LineWidth', 1.4)
    grid on
    
end