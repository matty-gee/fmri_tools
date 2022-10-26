function cond_struct = spm_cond_struct(name, onsets, durations, trials)

    %------------------------------------------------------------------------------
    % Creates a spm condition structure 
    %
    % Arguments
    % ---------
    % name : str
    % onsets : numeric array
    % durations : numeric array
    % TODO: add kwargs
    % [By Matthew Schafer, github: @matty-gee; 2020ish] 
    %------------------------------------------------------------------------------

    cond_struct = struct('name', name, 'onset', onsets, 'duration', durations,...
                         'trials', trials,...                         
                         'pmod', struct('name', {}, 'param', {}, 'poly', {}));