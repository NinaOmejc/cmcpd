function [data_h, data_p] = correct_for_uneven_size(data_h, data_p)
nP = size(data_p, 1); nH = size(data_h, 1);
if nH < nP 
    data_h = [data_h; NaN(nP-nH, 1)];
else
    data_p = [data_p; NaN(nH-nP, 1)];    
end
end