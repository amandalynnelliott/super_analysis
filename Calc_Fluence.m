function [Fluence, Fluence_norm] = Calc_Fluence(EF_arr,dt)

    Constants;

    % Calculate fluence by integrating out time, Fluence(r,z) -----------------
    Fluence = (1/2) * const.cl * const.eps0 * trapz(abs(EF_arr).^2, 2) * dt;
    Fluence = squeeze(Fluence);
    
    % Calculate fluence normalized in z ---------------------------------------
    temp_u = size(Fluence);
    Fluence_norm = zeros(temp_u);
    for i = 1:temp_u(1)                 % loop in z
        temp_max = max(Fluence,[],1);   % get max in that z
        Fluence_norm(i,:) = Fluence(i,:) ./ temp_max; 
    end

    clearvars temp_u temp_max
end