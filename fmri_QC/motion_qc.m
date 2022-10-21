clear

% navigate to correct directory, with sub/func structure
rpDir = uigetdir();
n = neuroelf;
rps = n.findfiles(rpDir, 'rp*.txt');

for f = 1:length(rps)

    driftT(f, :) = rmsDriftP20(rps{f});
    diffT(f, :) = rmsDiffP20(rps{f});
    
end

diffT = cell2table(diffT, 'VariableNames',{'SubID', 'RunID', 'MaxDiff', 'MeanDiff'});
writetable(diffT, 'totalstats_diff.txt');

driftT = cell2table(driftT, 'VariableNames',{'SubID', 'RunID', 'MaxDrift', 'MeanDrift'});
writetable(driftT, 'totalstats_drift.txt');