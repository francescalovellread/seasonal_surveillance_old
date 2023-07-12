function theta = symptoms(t)

    precip = [46.7 49 45.3 33.9 22.0 13.1 10.8 16.5 42.7 61.6 70.1 58.0];
    precip=precip-min(precip);
    precip=precip/max(precip);
    dryness = 1-precip;
    
    month = ceil(mod(t,365)*12/365);
    if month==0; month=12; end
    
    theta = dryness(month);

end

% function theta = symptoms(t)
%     
%     theta = 1;
% 
% end