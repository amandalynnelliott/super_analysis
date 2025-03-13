function [Q_greater, Q_lesser, temp, g_total] = Calc_SCG_Energy(omg0_1, omg0_2, ...
    omg_arr, spectrum)
    
    [~,i_1] = min(abs(omg_arr - omg0_1));
    [~,i_2] = min(abs(omg_arr - omg0_2));
    
    delta = 2 * (omg0_2 - omg0_1);
%     delta_index = ceil(delta / omg_arr); 
    d_omg = omg_arr(2) - omg_arr(1);
    delta_index = delta/d_omg;
    
    g_total = i_2 + delta_index;
    temp = i_1 - delta_index;
    
    Q_lesser = cumtrapz(spectrum(:,1:temp),2);
    Q_greater = cumtrapz(spectrum(:,g_total:end),2);
    Q_total = trapz(spectrum(:,1:end),2);
    
    u = size(spectrum);
    
    for iz=1:u(1)
        Q_lesser(iz,:) = Q_lesser(iz,:)/Q_total(iz);
        Q_greater(iz,:) = Q_greater(iz,:)/Q_total(iz);
    end

    clearvars iz i_1 i_2 delta delta_index u Q_total;
end