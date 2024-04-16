

rps = neuroelf.findfiles(pwd, 'rp*.txt');
rmsstats = cell(length(rps), 4);

for i = 1:length(rps)
    subid = regexp(rps{i}, '\<per\d\d\d\d', 'match');
    runid = regexp(rps{i}, '\<r\d\d', 'match');
    
    fname = strcat(subid, '_', runid, '_rms_diff_2019.txt');
    fname = char(fname);
    
    rparams = load(rps{i});
    rparams(:, 4:6) = 50 .* rparams(:, 4:6);%have to convert radians to millimeters
    rpdiff = zeros(length(rparams),6);
    
    for j = 2:length(rparams)
        diffcalc = rparams(j,:) - rparams(j-1,:);
        rpdiff(j,:) = diffcalc;
        rpdiff(1,:) = rparams(1,:);
    end
    
    rpdiffsq = rpdiff.^2;
    rpdifsqsum = sum(rpdiffsq,2);
    rpdifsumsqrt = rpdifsqsum.^.5;
    
    maxvalue = max(rpdifsumsqrt);
    meanvalue = mean(rpdifsumsqrt);
    maxvalue = num2cell(maxvalue);
    meanvalue = num2cell(meanvalue);
    
    rmsstats(i,1) = subid;
    rmsstats(i,2) = runid;
    rmsstats(i,3) = maxvalue;
    rmsstats(i,4) = meanvalue;
    
    save(fname, 'rpdifsumsqrt', '-ascii');
    
end

T = cell2table(rmsstats,'VariableNames',{'SubID', 'RunID', 'Max', 'Mean'});
writetable(T, 'totalstats_diff.txt');