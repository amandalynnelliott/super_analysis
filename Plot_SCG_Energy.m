% Calculate integrals -------------------------------------------------
[output.E_spectral.Q_greater, output.E_spectral.Q_lesser, ...
    temp, g_total] = ...
    Calc_SCG_Energy(pulse{1}.omg0,pulse{2}.omg0, grid.omg, ...
    shifted_Spectrum);

% Plot ----------------------------------------------------------------
figure;
imagesc(Spectrum.axes(1).data(g_total:end),Spectrum.axes(2).data, ...
    output.E_spectral.Q_greater); axis xy;
xlabel("$\omega [\omega_0]$"); ylabel("z [cm]"); 
title("$Q_{>}(\omega,z)$"); colorbar; 
filename = plots_folder + sim_title + "_Q_greater.jpg";
saveas(gcf,filename);
clearvars filename;

figure;
imagesc(Spectrum.axes(1).data(1:temp),Spectrum.axes(2).data, ...
    output.E_spectral.Q_lesser(:,1:temp)); axis xy;
xlabel("$\omega [\omega_0]$"); ylabel("z [cm]"); 
title("$Q_{<}(\omega,z)$"); colorbar; 
set(gca, 'ColorScale', 'log');
filename = plots_folder + sim_title + "_Q_lesser.jpg";
saveas(gcf,filename);
clearvars filename;

clearvars temp g_total;