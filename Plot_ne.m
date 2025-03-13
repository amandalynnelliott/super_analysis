if json_arr.laser.ionization
    % Electron density on axis --------------------------------------------
    neOA = SupData(data_folder + 'neOA.h5');
    
    figure; imagesc(grid.t, grid.z, abs(neOA.data).'); 
    axis xy;
    colorbar; colormap(hot); 
%     set(gca, 'ColorScale', 'log');
    xlabel("t [s]"); ylabel("z [cm]"); title("Electron Density on Axis")
    filename = plots_folder + sim_title + "_neOA.jpg";
    saveas(gcf,filename);
    clearvars filename;

    max_ne_init = max(abs(neOA.data(:,1)), [], 'all');
    disp(['max ne initial: ' num2str(max_ne_init)]);
    clearvars max_ne_init;
    
    % Electron density at each z step -------------------------------------
    if plot_arr.neF_full
        neF = SupData(data_folder + 'neF.h5');
        
        mkdir(plots_folder + "neF/")
        figure; 
        
        size_neF = size(neF.data);
        
        for iz = 1:size_neF(3)-1
            imagesc(neF.axes(2).data, neF.axes(1).data, ...
                squeeze(neF.data(:,:,iz))); axis xy; hold off;
            colorbar; colormap(hot); set(gca, 'ColorScale', 'log');
            a = colorbar;
            a.Label.String = "Electric Density";
            xlabel(neF.axes(2).name + neF.axes(2).units); 
            ylabel(neF.axes(1).name + neF.axes(1).units);
            title("z = " + num2str(neF.axes(3).data(iz)) + ...
                neF.axes(3).units);
            filename = plots_folder + "neF/" + sim_title + "_" + ...
                num2str(iz) + "_neF.jpg";
            saveas(gcf,filename);
            clearvars filename;
        end
        
        clearvars size_neF iz;
    end
end
