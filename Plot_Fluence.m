if exist('EF', 'var') ~= 1
    EF = SupData(data_folder + 'EF-re.h5',true);
end


[Fluence, Fluence_norm] = Calc_Fluence(EF.data, grid.dt);

% Plot fluence in r vs z (non-normalized) ---------------------------------
figure;
imagesc(grid.z,grid.r(1:length(grid.r)/6), Fluence(1:length(grid.r)/6,:)); 
axis xy;
xlabel(EF.axes(3).name + EF.axes(3).units); 
ylabel(EF.axes(1).name + EF.axes(1).units); title("Fluence [$J/m^2$]"); 
colorbar; 
filename = plots_folder + sim_title + "_fluence.jpg";
saveas(gcf,filename);
clearvars filename;


% Plot fluence in r vs z (normalized) ---------------------------------
figure;
imagesc(grid.z,grid.r(1:length(grid.r)/6), ...
    Fluence_norm(1:length(grid.r)/6,:)); 
axis xy;
xlabel(EF.axes(3).name + EF.axes(3).units); 
ylabel(EF.axes(1).name + EF.axes(1).units); 
title("Fluence Normalized"); colorbar; 
filename = plots_folder + sim_title + "_fluence_norm.jpg";
saveas(gcf,filename);
clearvars filename Fluence Fluence_norm;
