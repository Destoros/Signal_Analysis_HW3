function symbol = get_parameters(xn, fS)

%idea compute 2 spectrogram, one with a good time resolution, another one
%with a good frequency resoltuion

%FIRST: GOOD TIME RESOLUTION
window_length = 25;

[s, f, t] = spectrogram(xn, window_length, 0, 250, fS); %apply 0-padding
                                                        %nfft = 250 looked
                                                        %reasonable for the
                                                        %plot

s_dB = 20*log10(abs(s)); %calucalte signal level in dB
    
%SEEN FROM MESH PLOT:
%dB of signal is 22dB


%only look on time resolution
%find the places whre the signal level is higher than 20dB
bool_time = s_dB > 12; 

bool_time = logical(sum(bool_time,1));%add up all rows (remove frequency axis)



index = [];
flag = 0; 
for ii = 1:length(bool_time)
    %find the start and end indizes of a signal
    %at the start bool_time is 0 for sure
    
    if bool_time(ii) == 1 && flag == 0
        index = [index, ii]; %save current index
        flag = 1; %set flag to 1, so this loop searches for the next 0
    end
    
    if bool_time(ii) == 0 && flag == 1
        index = [index, ii]; 
        flag = 0;
    end
    
end


%at this time instances a signal changes
time_change_points = t(index); %we assume that the signal is 0 at the beginning (defined by flag = 0 before loop)

t_end = (length(xn) - 1000)/fS + window_length/2/fS; %time for t1 of last symbol
time_change_points = [time_change_points(:); t_end].';

t_diff = diff(time_change_points); %get the difference between the instances -> results t0 and t1



samples = t_diff*fS; %calculated back into samples

%round to nearest xx500 number
samples = round(samples*2/1000)*1000/2;




%SECOND: GOOD FREQUENCY RESOLUTION

window_length = 500; %window_length = 1000 has a too bad time resolution

[s, f, t] = spectrogram(xn, window_length, 0, 4*window_length, fS); %apply 0-padding

[T, ~] = meshgrid(t, f); %convert to mesh data

s_dB = 20*log10(abs(s));


s_threshold = 20; %in dB; seen from plot; results in several different frequencies for one symbol, need to be considered

bool_freq = s_dB > s_threshold;

%split into different symbols
n_symbols = length(samples)/2; %every symbol gives two sample numbers

frequency = zeros(n_symbols,1);

for ii = 1:n_symbols
    
    curr_time_boundaries = time_change_points([2*ii-1,2*ii]);
    
    T_slice = T >= curr_time_boundaries(1) & T <= curr_time_boundaries(2);
    
    F_slice = bool_freq & T_slice;
    
    f_slice = logical(sum(F_slice,2)); %without logical it could be used as a weighted sum
    
    frequency(ii) = mean(f(f_slice));
    
end

frequency = round(frequency,1);

for ii = 1:n_symbols
    symbol(ii).s0 = samples(2*ii-1);
    symbol(ii).s1 = samples(2*ii);
    symbol(ii).freq = frequency(ii);
end





end