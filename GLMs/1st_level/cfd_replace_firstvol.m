function cfd_num = cfd_replace_firstvol(cfd, method)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Replace the first volume in fmriprep confound by some method
    %
    % Arguments
    % ---------
    % cfd : cell array of strings
    % method: string
    %   'mean' or 'zero'
    %
    % [By Matthew Schafer, github: @matty-gee; 2020ish] 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s

    % 

    cfd_num = str2num(cfd(2:end,:)); % it's strings
    if strcmp(method, 'mean')
        cfd_num = [mean(cfd_num); cfd_num];
    elseif strcmp(method, 'zero')
        cfd_num = [0; cfd_num];
    end
    
end

