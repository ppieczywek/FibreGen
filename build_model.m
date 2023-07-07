%%%%%%%%%%%%%%%%%%%%%% define base model parametes  %%%%%%%%%%%%%%%%%%%%%%%
fibre_radius = 25;          %nm
segment_length = 25;        %nm
segments_number = 100;       %number of segments

model_width = 3000;         %nm
model_thickness = 1000;     %nm
model_length = 7000;        %nm
network_density = 0.4;      %percentage

%%%%%%%%%%%%%%%%%%%%%%%% generate fibre network  %%%%%%%%%%%%%%%%%%%%%%%%%%
[fibres] = generate_fibre_network(fibre_radius, segment_length, segments_number, ...
                                  0, 2*pi, 0.02, ...
                                  0, 0.01, 0.01, ...
                                  model_width, model_thickness, model_length, network_density);

%%%%%%%%%%%%%%%%%%%% define base units for simulations %%%%%%%%%%%%%%%%%%%%
% Using reduced units enables to avoid excessively small or big numbers
% being used during calculations.
kb = 1.380649 * 10^(-23);   
T  = 298;                   %K
base_energy = kb*T*100000;
base_length = 50;           %nm
base_volume = pi * 25^2 * (1.25*25); 

%%%%%%%%%%%%%%%%%%% define fibre mechanical properties %%%%%%%%%%%%%%%%%%%%
young_modulus = 125 * 10^9; %GPa
segment_volume = (pi*fibre_radius^2)*segment_length;

fibre_radius_reduced = fibre_radius / base_length;
segment_length_reduced = segment_length / base_length;
segment_volume_reduced = (pi*fibre_radius_reduced^2)*segment_length_reduced; 
segment_mass_reduced = 1.0; 
    
young_modulus_reduced = (young_modulus * ((base_length*10^(-9))^3) / base_energy);
linear_stiffness  = young_modulus_reduced * (pi*fibre_radius_reduced ^2) / segment_length_reduced;
bend_stiffness = young_modulus_reduced * ((pi*(fibre_radius_reduced *2)^4)/64) / (2 * segment_length_reduced);


% In the following lines of code, a finished model of the cellulose network 
% is created, which can be used in simulations. The model is defined to 
% reflect uniaxial tensile test conditions. The model is saved to a text
% file, in the format used by the FibreNet molecular dynamics simulation 
% software (https://github.com/ppieczywek/FibreNet). 
fibres = fibres';
bead_positions = cell2mat(fibres);
model_max_coord = max(bead_positions);
model_min_coord = min(bead_positions);

figure;
if ~isempty(fibres)
    for ii=1:10:length(fibres)
        fibre = fibres{ii};
        plot3(fibre(:,1), fibre(:,2), fibre(:,3), 'b-');
        hold on;
    end
end
axis('equal');


% In order to imitate presence of clamps holding oposite edges of the model
% first, beads from regions near the edges will be recognized and indicated
% by individual bead type definitions. In this case all beads closer than
% 500 nm  to one of the edges will be defined as "clamps".
clamp_1_coord = (model_min_coord(3) + 500) / base_length;
clamp_2_coord = (model_max_coord(3) - 500) / base_length;

% Model will be shifthed along X, Y and Z axis by 1000 nm, to give beads
% enough space from the walls of simulation box.
model_shift_x = 1000 / base_length;
model_shift_y = 1000 / base_length;
model_shift_z = 1000 / base_length;

%%%%%%%%%%%%%%%%%%%%%%%%% build model file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = waitbar(0.0,'building model file...');
fid = fopen('test_model.txt', 'w', 'n','UTF-8');
bead_id = 0;
fibre_id = 0;
bond_id = 0;
bead_group_id = 0;
bond_group_id = 0;
angle_group_id = 0;

for ii=1:1:length(fibres)

        fibre = fibres{ii} / base_length;
        fibre_id = ii;
        
        %%%%%%%%%%%%%%%%%%%%%% defines fibre beads %%%%%%%%%%%%%%%%%%%%%%%%
        for pp=1:1:size(fibre,1)
            if (fibre(pp,3) < clamp_1_coord)
                fibre(pp,1) = fibre(pp,1) + model_shift_x;
                fibre(pp,2) = fibre(pp,2) + model_shift_y;
                fibre(pp,3) = fibre(pp,3) + model_shift_z;
                fprintf(fid, 'BEAD\t%s\t%i\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\r\n','C3', ...
                    fibre_id, bead_group_id, fibre_radius_reduced, segment_mass_reduced, ...
                    fibre(pp,1), fibre(pp,2), fibre(pp,3), ...
                    0.0, 0.0, 0.0, ...
                    0.0, 0.0, 0.0, ...
                    fibre(pp,1), fibre(pp,2), fibre(pp,3));
                    
            elseif (fibre(pp,3) > (clamp_2_coord)) 
                fibre(pp,1) = fibre(pp,1) + model_shift_x;
                fibre(pp,2) = fibre(pp,2) + model_shift_y;
                fibre(pp,3) = fibre(pp,3) + model_shift_z;
                fprintf(fid, 'BEAD\t%s\t%i\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\r\n','C2', ...
                    fibre_id, bead_group_id, fibre_radius_reduced, segment_mass_reduced, ...
                    fibre(pp,1), fibre(pp,2), fibre(pp,3), ...
                    0.0, 0.0, 0.0, ...
                    0.0, 0.0, 0.0, ...
                    fibre(pp,1), fibre(pp,2), fibre(pp,3));
            else
                fibre(pp,1) = fibre(pp,1) + model_shift_x;
                fibre(pp,2) = fibre(pp,2) + model_shift_y;
                fibre(pp,3) = fibre(pp,3) + model_shift_z;
                fprintf(fid, 'BEAD\t%s\t%i\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\r\n','C1', ...
                    fibre_id, bead_group_id, fibre_radius_reduced, segment_mass_reduced, ...
                    fibre(pp,1), fibre(pp,2), fibre(pp,3), ...
                    0.0, 0.0, 0.0, ...
                    0.0, 0.0, 0.0, ...
                    fibre(pp,1), fibre(pp,2), fibre(pp,3));
            end
            bead_id = bead_id + 1;
        end
    
        %%%%%%%%%%%%%%%%%%%%%% defines bond potentials %%%%%%%%%%%%%%%%%%%%
        for pp = bead_id-size(fibre,1):1:bead_id-2 
            fprintf(fid, 'SPRING\t%s\t%i\t%i\t%i\t%.2f\t%.2f\r\n', ...
                'HR',bond_group_id, pp, pp+1, segment_length_reduced, linear_stiffness);            
        end
        
        %%%%%%%%%%% calculates rest angles for angle potentials %%%%%%%%%%%
        angles = [];
        for pp=2:1:size(fibre,1)-1

            ba = fibre(pp-1,:) - fibre(pp,:);
            bc = fibre(pp+1,:) - fibre(pp,:);
            l_ba = sqrt(sum(ba.^2));
            l_bc = sqrt(sum(bc.^2));
            dotp = (ba(1)*bc(1) + ba(2)*bc(2) + ba(3)*bc(3)) / (l_ba * l_bc);
            if dotp < -1.0 
               dotp = -0.999999;
            end
            if dotp >  1.0
               dotp =  0.999999;
            end
            angles = [angles; acos(dotp)];
        end
        
        %%%%%%%%%%%%%%%%%%%%% defines angle potentials %%%%%%%%%%%%%%%%%%%%
        angle_id = 1;
        for pp = bead_id-size(fibre,1):1:bead_id-3 
            fprintf(fid, 'ANGLE\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%.2f\t%.2f\r\n', ...
                'HR', angle_group_id, pp, pp+1, pp+2, bond_id, bond_id+1, angles(angle_id), bend_stiffness); 
            angle_id = angle_id + 1;
            bond_id = bond_id + 1;
        end
        bond_id = bond_id + 1;
        waitbar(ii/length(fibres), h, 'building model file...');
end
delete(h);
fclose(fid);
