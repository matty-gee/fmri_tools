function convert_spm_stat(conversion, infile, dof)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This script uses a template .mat batch script object to
    % convert images to a different statistical rep.
    % (Requires matlab stats toolbox)
    %
    %  Args:
    %  conversion -- one of 'TtoZ', 'ZtoT', '-log10PtoZ', 'Zto-log10P', 'PtoZ', 'ZtoP'
    %  infile -- input file stem (may include full path)
    %  outfile -- output file stem (may include full pasth)
    %  dof -- degrees of freedom
    %
    % Created by:           Josh Brown 
    % Modification date:    Aug. 3, 2007
    % Modified: 8/21/2009 Adam Krawitz - Added '-log10PtoZ' and 'Zto-log10P'
    % Modified: 2/10/2010 Adam Krawitz - Added 'PtoZ' and 'ZtoP'
    % Modifed: Ms - take nifties instead of .img files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmp(conversion,'TtoZ')
        expval = ['norminv(tcdf(i1,' num2str(dof) '),0,1)'];
    elseif strcmp(conversion,'ZtoT')
        expval = ['tinv(normcdf(i1,0,1),' num2str(dof) ')'];
    elseif strcmp(conversion,'-log10PtoZ')
        expval = 'norminv(1-10.^(-i1),0,1)';
    elseif strcmp(conversion,'Zto-log10P')
        expval = '-log10(1-normcdf(i1,0,1))';
    elseif strcmp(conversion,'PtoZ')
        expval = 'norminv(1-i1,0,1)';
    elseif strcmp(conversion,'ZtoP')
        expval = '1-normcdf(i1,0,1)';
    else
        disp(['Conversion "' conversion '" unrecognized']);
        return;
    end

    outfile = [erase(infile, '.nii') '_' conversion '.nii'];
    if ~isfile(outfile)
        jobs{1}.util{1}.imcalc.input{1} = [infile ',1'];
        jobs{1}.util{1}.imcalc.output = outfile;
        jobs{1}.util{1}.imcalc.expression = expval;
        spm_jobman('run', jobs);
        disp(['Conversion ' conversion ' complete.']);       
    else
        disp('File already exists.');
    end