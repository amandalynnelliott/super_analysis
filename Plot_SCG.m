
close all;

% Find the 2nd and 3rd harmonics
harm2 = 2 * pulse{1}.omg0;

try 
    harm3 = 3 * pulse{2}.omg0;
catch 
    harm3 = 3 * pulse{1}.omg0;
end

% Get the midpoint between them
harm_mid = (harm2 + harm3)/2;

% Find the nearest index
[~,harm_mid_index] = min(abs(harm_mid - grid.omg));

% Get index of the focal position, zf.
zf = json_arr.laser.pulses(1).zf * 100; % m -> cm
[~,zf_index] = min(abs(zf - Spectrum.axes(2).data));

% Calculate the SCG
SCG_Spectrum = Calc_Supercontinuum_Spectra(shifted_Spectrum,zf_index);

max_SCG_spectrum = max(max(SCG_Spectrum));

figure; 
imagesc(Spectrum.axes(1).data(i_zero:harm_mid_index),Spectrum.axes(2).data,...
    SCG_Spectrum(:, i_zero:harm_mid_index)); 
axis xy; 
colorbar; 
set(gca, 'ColorScale','log'); 
% set(gca, 'CLim', [orders*max_SCG_spectrum,max_SCG_spectrum]); % Orders of magnitude
colormap(seismic);
xlabel("$\omega [\omega_0]$"); 
ylabel("$z$ [cm]"); title('SCG Spectrum Method 1b')
filename = plots_folder + sim_title + "_Spectrum_SCG1-modified.jpg";
saveas(gcf,filename);
% clearvars filename harm2 harm3 harm_mid_index harm_mid i_zero zf zf_index;