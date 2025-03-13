% outer_folder = "h-runs_sub/";
% load(outer_folder + 'energy_constraint')


target_U = 1.996e-4;
target_obj = [init_U_arr_super, init_U_arr_sub];

for i=1:length(target_obj)
    for j=1:length(target_obj(i).energy)
        percent_diff = abs(target_obj(i).energy(j) - target_U)/target_U; 

        target_obj(i).percent_diff(j) = percent_diff;

        if percent_diff >= 0.05
            disp(target_obj(i).sim(j))
            disp(target_obj(i).percent_diff(j))
        end
    end
end