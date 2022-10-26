function cfd_json = fmriprep_cfd_json_read(cfd_json_fname)
    
    %------------------------------------------------------------------------------
    % Opens fmriprep confound json and makes readable in matlab
    %
    % Arguments
    % ---------
    % cfd_json_fname : str
    %   path to find fmriprep confound json
    %
    % [By Matthew Schafer, github: @matty-gee; 2020ish] 
    %------------------------------------------------------------------------------s

    fid = fopen(cfd_json_fname); 
    raw = fread(fid, inf); 
    str = char(raw'); 
    fclose(fid);
    cfd_json = jsondecode(str); 