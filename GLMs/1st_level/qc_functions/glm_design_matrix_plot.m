% design matrix plotting
events   = {glm.cond(2).onset};
n_vols   = 1570; % number of overall volumes
n_events = length(events{1});   % number of events


% % Create and display a table of these onset times:
% t = table(onsets{1}(:, 1), onsets{1}(:, 2), onsets{2}(:, 1), onsets{2}(:, 2), 'VariableNames', {'Evt1_Time' 'Evt1_Dur' 'Evt2_Time' 'Evt2_Dur'});
% disp(t)
% 
% clear X bfname

% Build three models: Convolve with three basis sets (need SPM on path)
% 1: canonical, 2: 3-parameter, 3: FIR

bfname{1} = 'Canonical HRF';
X{1} = onsets2fmridesign(events, 1, n_vols, spm_hrf(1));

bfname{2} = 'hrf (with time and dispersion derivatives)';
X{2} = onsets2fmridesign(events, 1, n_vols, bfname{2});

bfname{3} = 'Finite Impulse Response';
X{3} = onsets2fmridesign(events, 1, n_vols, bfname{3});

create_figure('X_3basis sets', 2, 3);

for i = 1:3

    subplot(2, 3, i)
    h = plot_matrix_cols(zscore(X{i}(:, 1:end-1)), 'vertical');
    set(gca, 'XColor', 'w', 'YTick', [0:20:100]);
    axis tight
    title(bfname{i})
    ylabel('Time')

    subplot(2, 3, 3+i)
    imagesc(X{i}(:, 1:end-1))
    set(gca, 'YDir', 'Reverse', 'XTickLabel', []);
    ylabel('Time')
    axis tight

end

drawnow, snapnow