% Initial E-field ---------------------------------------------------------
Einit = SupData(data_folder + 'Einit-re.h5',true); 

grid.r = Einit.axes(1).data;
grid.dr = Einit.axes(1).dax;
grid.t = Einit.axes(2).data;
grid.dt = Einit.axes(2).dax;

figure; imagesc(grid.t, grid.r, abs(Einit.data)); 
axis xy; 
colorbar; colormap(parula);
xlabel("t [s]"); ylabel("r [m]"); title("Initial Field")
filename = plots_folder + sim_title + "_Einit.jpg";
saveas(gcf,filename);
clearvars filename;

max_Einit = max(abs(Einit.data), [], 'all');
disp(['max Einit: ' num2str(max_Einit)]);
clearvars max_Einit;

% Final E-field -----------------------------------------------------------
try 
    Efinal = SupData(data_folder + 'Efinal-re.h5',true); 
catch 
    warning('Efinal is undefined. E-field may have gone to NaNs.')
end

if exist('Efinal', 'var')
    figure; imagesc(grid.t, grid.r, abs(Efinal.data)); 
    axis xy; 
    colorbar; colormap(parula);
    xlabel("t [s]"); ylabel("r [m]"); title("Final Field")
    filename = plots_folder + sim_title + "_Efinal.jpg";
    saveas(gcf,filename);
    clearvars filename;
end

% E-field on axis ---------------------------------------------------------
EOA = SupData(data_folder + 'EOA-re.h5',true);

grid.z = EOA.axes(2).data * 100;
grid.dz = EOA.axes(2).dax * 100;

figure; imagesc(grid.t, grid.z, abs(EOA.data).'); axis xy;
colorbar; colormap(parula);
xlabel("t [s]"); ylabel("z [cm]"); title("E-Field on Axis")
filename = plots_folder + sim_title + "_EOA.jpg";
saveas(gcf,filename);
clearvars filename;

% E-field on axis lineout at pulse center ---------------------------------
zf = json_arr.laser.pulses(1).zf;
[~,zindex] = min(abs(zf - EOA.axes(2).data));
[~, tindex] = max(abs(EOA.data(:,zindex)));
EOA_slice = abs(EOA.data(tindex,:));

figure; 
plot(grid.z, EOA_slice); 
xlabel("z [cm]"); ylabel("[V/m]"); 
title("E-Field on Axis at t=0"); 
filename = plots_folder + sim_title + "_EOA-lineout.jpg";
saveas(gcf,filename);
clearvars filename zf zindex tindex EOA_slice;

% E-field at each z step --------------------------------------------------
if plot_arr.Efield_full
    EF = SupData(data_folder + 'EF-re.h5',true);
    
    mkdir(plots_folder + "EF/")
    figure; 
    
    clearvars iz;
    size_EF = size(EF.data);

    if plot_arr.GIF
        max_EF = max(abs(EF.data), [], 'all');
        min_EF = min(abs(EF.data), [], 'all');
    end
    
    
    for iz = 1:size_EF(3)-1
        imagesc(EF.axes(2).data, EF.axes(1).data, ...
            abs(squeeze(EF.data(:,:,iz)))); axis xy; hold off;
        colorbar; colormap(parula); 
        if plot_arr.GIF
             set(gca, 'CLim', [min_EF, max_EF]);
        else
            set(gca, 'ColorScale', 'log');
        end

        a = colorbar;
        a.Label.String = "Electric Field";
        xlabel(EF.axes(2).name + EF.axes(2).units); 
        ylabel(EF.axes(1).name + EF.axes(1).units);
        title("z = " + num2str(EF.axes(3).data(iz)) + EF.axes(3).units);
        filename = plots_folder + "EF/" + sim_title +...
            "_" + num2str(iz) + "_EF.jpg";
        
        saveas(gcf,filename);
        clearvars filename;
    end
    
    clearvars size_EF iz max_EF min_EF;
end