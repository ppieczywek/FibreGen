function save_model_to_xyz(fibers, file_name)
%
%   Function takes coordinates of fibre network and saves them into XYZ
%   file.
%  
%
%   INPUT PARAMETERS 
%   fibers    - cell array, where each cell holds Nx3 matrix of coordinates
%               for each subsequent bead/segment of fibre;
%   file_name - name of the output XYZ file
 
    h = waitbar(0.0,'Saving data to XYZ...');
    fid = fopen( file_name,'w', 'n','UTF-8');
    patricles_num = sum(cellfun('length',fibers));    
    fprintf(fid, '%i\r\n', patricles_num );
    fprintf(fid, 'Model file generated with MATLAB script \r\n');

    for ii=1:1:length(fibers)
        for jj=1:1:length(fibers{ii})
            fprintf(fid, 'C\t%f\t%f\t%f\r\n', fibers{ii}(jj,1), fibers{ii}(jj,2), fibers{ii}(jj,3) );
        end
        waitbar(ii/length(fibers), h, 'Saving data to XYZ...');
    end
    fclose(fid);    
    delete(h);
   
end

