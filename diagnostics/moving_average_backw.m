% calculate the moving average of a 1D array with a specified width, only
% in the backwards direction. 


function out = moving_average_backw(input, width)
    out = NaN(size(input));
    for i = 1+width:length(input)
      
        temp = mean(input(i-width:i)); % calc it backwards
        out(i) = temp;
    end
end


