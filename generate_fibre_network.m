function [ fibre_coordinates ] = generate_fibre_network(radius, segment_length, ...
               segment_number, theta_start, theta_range, theta_dev, ...
               fi_start, fi_range, fi_dev, x_limits, y_limits, z_limits, target_density)
    
%   Function generates the coordinates of individual fiber segments in 3D 
%   space. The model fibers consist of a finite number of segments of equal 
%   length. The function allows you to control the fiber radius, length, 
%   position, alignment direction, curvature and packing density of the 
%   fibers in 3D space.
%
%   INPUT PARAMETERS 
%
%   radius              - positive floating point number describing the radius 
%                         of fibres
%   segment_length      - positive floating point number describing the
%                         distance between two consecutive beads
%   segment_number      - number of beads/segments per fibre
%   theta_start         - initial mean azimuthal angle of fibre;
%   theta_range         - the deviation of the azimuthal angle from the mean;  
%                         both theta_start and theta_range define the direction
%                         of the fibre;
%   fi_start            - inital polar angle of fibre;
%   fi_range            - the deviation of the polar angle from the mean;
%                         both theta_start and theta_range define the direction
%                         of the fibre;
%   theta_dev, fi_dev   - deviations in azimuthal and polar angles; defines the
%                         curvature of fibre; 
%   x_limits, 
%   y_limits, 
%   z_limits            - size of the network
%   target_density      - maximal fraction of volume occupied by fibres
%
%   Final output data is written to "fibre_ccordinats" cell array, where each 
%   cell holds Nx3 matrices of CG bead coordinates for each subsequent bead/segment
%   of the individual fibre.

    h = waitbar(0.0,'generating fibre network...');

    segment_volume = (pi*radius^2)*segment_length;

    sample_size = [x_limits y_limits z_limits]; 
    sample_volume =  sample_size(1)*sample_size(2)*sample_size(3);

    bias_x = ceil(max(sample_size)*0.2);
    bias_y = ceil(max(sample_size)*0.2);
    bias_z = ceil(max(sample_size)*0.2);

    model_margin_x = ceil(max(sample_size)*0.1);
    model_margin_y = ceil(max(sample_size)*0.1);
    model_margin_z = ceil(max(sample_size)*0.1);

    grid_min_x = -bias_x;
    grid_min_y = -bias_y;
    grid_min_z = -bias_z;

    grid_max_x = sample_size(1) + (2*model_margin_x) + bias_x;
    grid_max_y = sample_size(2) + (2*model_margin_y) + bias_y;
    grid_max_z = sample_size(3) + (2*model_margin_z) + bias_z;

    model_min_x = model_margin_x;
    model_min_y = model_margin_y;
    model_min_z = model_margin_z;

    model_max_x = sample_size(1) + model_margin_x;
    model_max_y = sample_size(2) + model_margin_y;
    model_max_z = sample_size(3) + model_margin_z;

    number_of_fibers = 1000000; 

    fibers = cell(1);
    inserted_fibers = 0;
    total_segments = 0;

    for qq=1:number_of_fibers

        posX = randi([grid_min_x, grid_max_x],1);
        posY = randi([grid_min_y, grid_max_y],1);
        posZ = randi([grid_min_z, grid_max_z],1);

        
        
        
        segments_num = segment_number;
        theta_start = theta_start + theta_range*(rand(1)-0.5);
        fi_start = fi_start + fi_range*(rand(1)-0.5);
        
        segments = [posX posY posZ];

        for ii=1:1:segments_num

                theta = theta_start + theta_dev*(rand(1)-0.5);
                theta_start = theta;

                fi = fi_start + fi_dev*(rand(1)-0.5);
                fi_start = fi;

                u(1) = cos(theta)*cos(fi);
                u(2) = cos(theta)*sin(fi);
                u(3) = sin(theta);

                x = posX + segment_length*u(1);
                y = posY + segment_length*u(2);
                z = posZ + segment_length*u(3);

                if x < grid_min_x || x > grid_max_x
                    break;
                elseif y < grid_min_y || y > grid_max_y
                    break;
                elseif z < grid_min_z || z > grid_max_z
                    break;
                else            
                    segments = [segments; [x y z]];
                    posX = x;
                    posY = y;
                    posZ = z;
                end

        end

        ex_x = ((segments(:,1) < model_min_x) | (segments(:,1) > model_max_x));
        ex_y = ((segments(:,2) < model_min_y) | (segments(:,2) > model_max_y));
        ex_z = ((segments(:,3) < model_min_z) | (segments(:,3) > model_max_z));
        segments(ex_x | ex_y | ex_z,:) = [];   

        if ~isempty(segments)
            max_segments_length = max(sqrt(sum(diff(segments).^2,2)));
            if size(segments,1) > 15 && max_segments_length < (1.1*segment_length)
                inserted_fibers = inserted_fibers + 1;
                segments(:,1) = segments(:,1) - model_min_x;
                segments(:,2) = segments(:,2) - model_min_y;
                segments(:,3) = segments(:,3) - model_min_z;
                fibers{inserted_fibers} = segments;
                total_segments = total_segments + size(segments,1) * segment_volume;
                if (total_segments/sample_volume) >= target_density 
                    break;
                end
                waitbar((total_segments/sample_volume)/target_density  ,h,'generating fibre network...');
            end
        end
    end

    delete(h)
    fibre_coordinates = fibers;
end
