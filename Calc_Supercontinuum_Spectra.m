function [SCG_Spectrum] = Calc_Supercontinuum_Spectra(Spectra,zf_index)
% Focus on the supercontinuum generation close around the initial
% bandwidth.
% However, the original spectra (z=0) is subtracted.
% Input is the shifted spectra from Streak.h5, with (omg, z) dimensions.
    u = size(Spectra);
    SCG_Spectrum = zeros(u);
    init_spectra = Spectra(1,:);  % Method 1 and 3
%     figure; plot(init_spectra);
%     init_spectra = mean(Spectra(1:zf_index,:)); % Method 2

% Methods 1 and 2
    for iz = 1:u(1)
        SCG_Spectrum(iz,:) = Spectra(iz,:) - init_spectra;
    end

    % ---------------------------------------------------------------------
    % SYMLOG---------------------------------------------------------------
    % ---------------------------------------------------------------------
    
        % Applying symlog function to deal with plotting negative values
    % in log scaling
%     C = 5; % scaling constant that determines the resolution of the data
            % around zero. The smallest order of magnitude on either side
            % of zero will be 10^ceil(C).
%     SCG_Spectrum_symlog = sign(SCG_Spectrum).*(log10(1 + ...
%         abs(SCG_Spectrum)./(10^C)));

%     spectra_lineout = Spectra(zf_index,:);
%     SCG_spectra_lineout = SCG_Spectrum(zf_index,:);
%     figure; plot(spectra_lineout); hold on; plot(SCG_spectra_lineout); hold off; title('Before symlog');
%     legend('spectra', 'SCG'); set(gca, 'YScale','log'); 
% 
%     spectra_lineout_symlog = sign(spectra_lineout).*(log10(1 + ...
%         abs(spectra_lineout)./(10^C)));
%     SCG_spectra_lineout_symlog = sign(SCG_spectra_lineout).*(log10(1 + ...
%         abs(SCG_spectra_lineout)./(10^C)));
%     figure; plot(spectra_lineout_symlog); hold on; plot(SCG_spectra_lineout_symlog); hold off; title('After symlog');
%     legend('spectra', 'SCG'); set(gca, 'YScale','log'); 

%     figure; plot(SCG_Spectrum(1,:))

    % ---------------------------------------------------------------------
    % end SYMLOG-----------------------------------------------------------
    % ---------------------------------------------------------------------

% Method 3

%     % Cast nonzero values to 1 and ~zero values to 0.
%     order = 1e-4;
%     init_Spectra_mod = init_spectra;
%     size(init_Spectra_mod)
%     init_Spectra_mod(init_Spectra_mod <= max(init_spectra)*order) =0;
%     init_Spectra_mod(init_Spectra_mod >  max(init_spectra)*order) =1;
% 
% %     figure; plot(init_Spectra_mod);
% 
%     % Mask the rest of the spectra
%     SCG_Spectrum = Spectra;
%     SCG_Spectrum(:, init_Spectra_mod == 1) = min(min(Spectra));
    
    clearvars u iz init_spectra
end