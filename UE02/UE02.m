close all
clear all
clc


id = 25; %personalized example
fS = 8192; %Hz


%(A) ----------------------------------------------------------------------
% Decoding using STFT

xn = signalgenerator2(id);
% plot(xn)

figure
    %s = spectrogram(x,window,noverlap,nfft,fs)
    spectrogram(xn, 500, 0, 500)    %window: length of window;
                                    %noverlap: number of overlap
                                    %nfft: n-point DFT
                                    %fs: sampling frequency 
    title(['Spectogramm of given signal for id =' num2str(id)])
    
                                
                               
                                
%I chose window length = 500  because its the gcd of the given sample numbers for the 2nd and 3rd
%encryption length (and also a divided of the 1000 pre and post appened
%zero samples)


%READ FROM SPECTROMGRAMM PLOT USING DATA CURSOR:
%DUE TO THE WAY SPECTOGRAM WORKS THE, ALL SAMPLES ARE DISPLAYED WITH AN 
%OFFSET OF 250, BUT THIS IS ONLY IMPORTANT FOR THE LAST t1
%First symbol: 
%   f = 0.12*fS/2 = 491 Hz                      -> ' '
%   t0 = 4250 - 750 = 3500                      -> 'z'
%   t1 = 14250 - 4250 = 10000                   -> 'w'

%Second symbol:
%   f = 0.288*fS/2 = 1176Hz                     -> 'e'
%   t0 = 1.525E4 - 1.425E4 = 1000               -> 'i'
%   t1 = 1.975E4 - 1.525E4 = 4500               -> ' '

%Third symbol:
%   f = 0.128*fS/2 = 524 Hz                     -> 'a'
%   t0 = 3.125E4 - 1.975E4 = 11500              -> 'c'
%   t1 = 3.375E4 - 3.125E4 = 2500               -> 'h'

%Fourth symbol:
%   f = 0.064*fS/2 = 262 Hz                     -> 't'
%   t0 = 3.825E4 - 3.375E4  = 4500              -> ' '
%   t1 = 4.375E4 - 3.825E4 - 1000 = 4500        -> ' '


%the 250 offset samples were considered in the last caluclation 



%(B) ----------------------------------------------------------------------
% Occupied bandwidth


%Time resolution is not as important as frequency resolution right now ->
%increase window length

window_length = 500; %Using window lengths which do not correspond to
% the possible delay times, the spectrogram delivers worse output 
%e.g. for window_length = 1000; the 2nd symbol which is very short, smears
%out a lot
%e.g. window_length = 5000; some symbols are almost not visible anymore
%e.g. window_length = 100; short window length -> bad frequency resolution
%window_length = 500 seems to be the best window.



[s, f, t] = spectrogram(xn, window_length, 0, 4*window_length, fS); %apply 0-padding
[T, F] = meshgrid(t, f); %convert to mesh data
s_dB = 20*log10(abs(s));



figure
    mesh(F,T,s_dB)
    xlabel('Frequency in Hz')
    ylabel('Time in s')
    zlabel('|X(jw)| in dB')
    title('Plot to find frequency boundaries')
    
    
% READ FROM PLOT
%Main Signal has ca. 42dB/rad/sample
%boundary = 42 - 20 = 22dB/rad/sample
highest_f = 1204; %Hz
lowest_f = 233.5; %Hz

bandwidth = highest_f - lowest_f; 
display(bandwidth);

%OTHER --------------------------------------------------------------------
% highest_f = 0.296*fS/2; % = 1212.4
% lowest_f = 0.056*fS/2; % = 229.4
%--------------------------------------------------------------------------



%(C) ----------------------------------------------------------------------
% generate table parameters 
string = 'abcdefghijklmnopqrstuvwxyz '; %vary through all possible characters
for ii = 1:length(string)
    xn = codingtest(string(ii));
%     spectrogram(xn, 500, 0, 500)
    symbol(ii) = get_parameters(xn, fS); %get both parameters for all characters
end

freq = [symbol.freq].';
s0 = [symbol.s0].'; %its in samples and not in time: t0 -> s0


%create table
col_names = {'Frequency', 'Length'};
row_names = cellstr(string.');
row_names{end} = '_';


T = table(freq, s0,'VariableNames', col_names, 'RowNames', row_names);
disp(T);
%The underscore represents the space character; RowNames need 
%non empty characters and a blank space counts as empty 


%save as excel file
T_save = table(string.', freq, s0,'VariableNames', ['character', col_names]);
writetable(T_save,'table.xlsx') %to easily include in presentation


%NOW GENERATE SIGNAL

soundstring1 = 'ghtqhtljtlhtxtnyaayhtxtnxht ht ht htaaa';
soundstring2 = 'xmmjmmymmktnkkm mmatnokmjmmymmommztnzkmkmm tnakm';


xn = own_encoder(soundstring1,fS, T_save); %Turbobier - Arbeitslos (https://www.youtube.com/watch?v=oRTQ1iK8fAY)
                                           %bzw Helene Fischer - Atemlos
% xn = own_encoder(soundstring2,fS, T_save); %Udo Jügrens - Ich war noch niemals in New York

[s, f, t] = spectrogram(xn, window_length, 0, 4*window_length, fS); %apply 0-padding
[T, F] = meshgrid(t, f); %convert to mesh data
s_dB = 20*log10(abs(s));



figure
    mesh(F,T,s_dB)
    xlabel('Frequency in Hz')
    ylabel('Time in s')
    zlabel('|X(jw)| in dB')
    title('Task 4: Plot to find frequency boundaries')      
    
    
%soundstrng 1:
% fmax = 1073;
% fmin = 467;

%soundstrng 2:
% fmax = 811;
% fmin = 364.5;






% (ADDITIONAL) ------------------------------------------------------------
%For comparison decode the given signal in Task 1 using the created function

xn = signalgenerator2(id);
string = own_decoder(xn, fS, T_save);
display(string)




a = 1; %just to set a breakpoint at the end of the file




