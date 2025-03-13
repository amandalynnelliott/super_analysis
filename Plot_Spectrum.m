Spectrum = SupData(data_folder + 'streak.h5');

grid.omg = Spectrum.axes(1).data;
grid.omg_norm = grid.omg / pulse{1}.omg0; % Normalize frequency
grid.domg = grid.omg(2) - grid.omg(1);

[~, i_zero] = min(abs(grid.omg/pulse{1}.omg0));

shifted_Spectrum = fftshift(Spectrum.data,1).';
max_spectrum = max(max(shifted_Spectrum));

if max_spectrum == inf || isnan(max_spectrum)
    max_spectrum = 1e3;
end

figure; 
imagesc(grid.omg_norm(i_zero:end),grid.z, shifted_Spectrum(:, i_zero:end)); 
axis xy; 
colorbar; 
set(gca, 'ColorScale','log'); 
set(gca, 'CLim', [orders*max_spectrum,max_spectrum]); % Orders of magnitude
colormap(parula);
xlabel("$\omega [\omega_0]$"); 
ylabel("$z$ [cm]"); title('Spectrum')
filename = plots_folder + sim_title + "_Spectrum.jpg";
saveas(gcf,filename);
clearvars filename max_spectrum;