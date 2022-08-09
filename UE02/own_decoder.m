function string = own_decoder(xn, fS, table)
%table entries should look like this (first row should be dismissed)
% 'character'     'frequency'     'Length/Pause'
%    a               523.3            5000
%   ...               ...              ...

symbol = get_parameters(xn, fS);
    
    
characters = table{:,1};
frequencies = table{:,2};
lengths = table{:,3};

string = [];
for ii = 1:length(symbol)
     %get the minmal squared distance
    [~, pos1] = min(abs(frequencies - symbol(ii).freq).^2);
    [~, pos2] = min(abs(lengths - symbol(ii).s0).^2);
    [~, pos3] = min(abs(lengths - symbol(ii).s1).^2);
    
    string = [string, characters([pos1, pos2, pos3]).'];
end


end