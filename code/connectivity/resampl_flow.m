% resampling of the t_series by moving average
%
% org_freq - sampling frequency of the original signal
%
% res_freq - sampling frequency of the resampled signal
%
% (t_series, 400, 10)

function [output] = resampl_flow(t_series, org_freq, res_freq)

dim = length(t_series);
    
T = org_freq / res_freq; 

output = zeros(floor(dim/T),1);
        
   for j = 1 : floor(dim/T)
        output(j) = mean( t_series( (j-1)*T+1 : (j-1)*T + T ) );
   end
    
   [dimx dimy]=size(output);

if dimx>dimy
    output = output';
end
% time = linspace(0, (floor(dim/T)-1) * 1/res_freq, floor(dim/T) ); 
   
