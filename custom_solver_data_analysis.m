% t_list: list of output time points, dimension: {number_of_intervals+1,1}
% conc_list: cell arrays of concentrations of output species, dimension: {number_of_intervals+1, number_of_species}
% conc_list: cell arrays of names of output species, dimension: {1, number_of_species}
% x: x coordinates of grid points, dimension: {ny, nx}
% y: y coordinates of grid points, dimension: {ny, nx}
% dx: the x direction length of grids, dimension: {ny, nx}
% dy: the y direction length of grid, dimension: {ny, nx}

function [t_list, conc_list, species_list, x, y, dx, dy] = custom_solver_data_analysis(number_of_intervals)
    eval('coordinates');
    x=vec2mat(eval('x_volume_global'),nx);
    y=vec2mat(eval('y_volume_global'),nx);
    dx=vec2mat(eval('dx_volume_global'),nx);
    dy=vec2mat(eval('dy_volume_global'),nx);
    
    species_list_file = fopen('output_species_list.txt','r');

    tline = fgetl(species_list_file);
    number_of_volume_species = str2double(tline);
    
    j=1;
    volume_species_list = cell(1,number_of_volume_species);
    while j<=number_of_volume_species && ischar(tline)
        tline=fgetl(species_list_file);
        volume_species_list{1,j} = tline;
        j=j+1;
    end
            
    if j~=number_of_volume_species+1
        disp('error, number of volume species does not agree in output_species_list.txt')
    end
    
    tline = fgetl(species_list_file);
    number_of_surface_species = str2double(tline);
    
    j=1;
    surface_species_list = cell(1,number_of_surface_species);
    while j<=number_of_surface_species && ischar(tline)
        tline=fgetl(species_list_file);
        surface_species_list{1,j} = tline;
        j=j+1;
    end
            
    if j~=number_of_surface_species+1
        disp('error, number of surface species does not agree in output_species_list.txt')
    end
    
    fclose(species_list_file);
   
    number_of_species = number_of_volume_species + number_of_surface_species;
    species_list = [volume_species_list, surface_species_list];

    t_list = cell(number_of_intervals+1,1);
    conc_list = cell(number_of_intervals+1,number_of_species);
    
    for i=1:number_of_intervals+1
        eval(['concentrations_' num2str(i-1)]);
        t_list{i} = t;
        for j=1:number_of_species
            conc_list{i,j} = vec2mat(eval(species_list{j}),nx);
        end
    end
