function [full_filter] = Calc_THz_filter(filter_gaussian, filter_half, freq_min, freq_max, grid_THz)

    switch filter_gaussian
        case '1' % Hard array, only 0's and 1's
            full_filter = (grid_THz >= freq_min) & (grid_THz <= freq_max);
    
        case '2' % Gaussian filter
            filter_ind = (grid_THz >= freq_min) & (grid_THz <= freq_max);
            first_ind = find(filter_ind,1);
            last_ind = find(filter_ind, 1, 'last');
            Nt_filt = last_ind - first_ind;
            
            [freq, filter] = gaussian_filter([freq_min,freq_max], (freq_min + freq_max)/2, 35,Nt_filt);
            figure; plot(freq,filter); xlabel('Freq (THz)'); ylabel('Filter Amplitude')
            
            temp_first = zeros(1, first_ind);
            temp_last = zeros(1, length(grid_THz) - last_ind);
            full_filter = [temp_first, filter, temp_last];
    
            clearvars Nt_filt temp_first temp_last;
    
            if filter_half  % Make it a hard cut-off at zero.
                middle_ind = ceil((last_ind + first_ind) /2);
                zero_ind = find(grid_THz == freq_min);

                full_filter(zero_ind+1:middle_ind) = 1;
               
                clearvars middle_ind zero_ind;
            end
    
            figure; plot(grid_THz,full_filter); xlabel('Freq (THz)'); ylabel('Filter Amplitude')
    
            clearvars filter_ind first_ind filter_gaussian filter_half;
    end

end