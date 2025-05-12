% Load participant data
participants = readtable('Participants.tsv', 'FileType', 'text', 'Delimiter', '\t');

% Initialize variables for averaging and MDS
modalities = {'haptic', 'visual'};
avgSimilarityMatrices = cell(1, 2);
stimuliSets = cell(1, 2);
validSubjectCounts = zeros(1, 2);
mdsConfigurations = cell(1, 2);
stressValues = zeros(1, 2);
stressByDimension = zeros(10, 2); % Stress for dimensions 1 to 10

for modIdx = 1:length(modalities)
    modality = modalities{modIdx};
    % Filter subjects for current modality
    modality_subjects = participants(strcmp(participants.Modality, modality), :);
    subjectIDs = modality_subjects.Participants;
    
    % Initialize variables for averaging
    allStimuli = {};
    sumSimilarityMatrix = [];
    nValidSubjects = 0;
    
    % Loop over each subject
    for sub = 1:length(subjectIDs)
        % Construct filename
        subjectID = subjectIDs{sub};
        filename = sprintf('%s_task-similarityRating_modality-%s.tsv', subjectID, [upper(modality(1)) modality(2:end)]);
        filename = fullfile(subjectID,filename)
        % Check if file exists
        if ~isfile(filename)
            fprintf('File %s not found, skipping...\n', filename);
            continue;
        end
        
        % Load data
        data = readtable(filename, 'FileType', 'text', 'Delimiter', '\t');
        
        % Extract unique stimuli
        stimuli = unique([data.FirstObj; data.SecondObj]);
        
        % Initialize similarity matrix
        nStimuli = length(stimuli);
        similarityMatrix = zeros(nStimuli, nStimuli);
        
        % Populate similarity matrix
        for i = 1:height(data)
            rowIdx = find(strcmp(stimuli, data.FirstObj{i}));
            colIdx = find(strcmp(stimuli, data.SecondObj{i}));
            similarityMatrix(rowIdx, colIdx) = data.ParticipantResp(i);
            similarityMatrix(colIdx, rowIdx) = data.ParticipantResp(i); % Symmetric matrix
        end
        
        % If first valid subject, set up stimuli and matrix
        if isempty(allStimuli)
            allStimuli = stimuli;
            sumSimilarityMatrix = zeros(nStimuli);
        end
        
        % Ensure stimuli are consistent
        if ~isequal(sort(allStimuli), sort(stimuli))
            fprintf('Inconsistent stimuli for %s, skipping...\n', subjectID);
            continue;
        end
        
        % Align similarity matrix to common stimuli order
        alignedMatrix = zeros(length(allStimuli));
        for i = 1:length(stimuli)
            for j = 1:length(stimuli)
                rowIdx = find(strcmp(allStimuli, stimuli{i}));
                colIdx = find(strcmp(allStimuli, stimuli{j}));
                alignedMatrix(rowIdx, colIdx) = similarityMatrix(i, j);
            end
        end
        
        % Accumulate for averaging
        sumSimilarityMatrix = sumSimilarityMatrix + alignedMatrix;
        nValidSubjects = nValidSubjects + 1;
    end
    
    % Store average matrix and stimuli
    if nValidSubjects > 0
        avgSimilarityMatrices{modIdx} = sumSimilarityMatrix / nValidSubjects;
        stimuliSets{modIdx} = allStimuli;
        validSubjectCounts(modIdx) = nValidSubjects;
        
        % Convert similarity to dissimilarity (1 = dissimilar, 7 = identical)
        dissimilarityMatrix = 8 - avgSimilarityMatrices{modIdx}; % Transform to dissimilarity
        dissimilarityMatrix(1:nStimuli+1:end) = 0; % Set diagonal to 0
        
        % Perform MDS for 2D configuration
        [Y, stress] = mdscale(dissimilarityMatrix, 2, 'Criterion', 'stress');
        mdsConfigurations{modIdx} = Y;
        stressValues(modIdx) = stress;
        fprintf('%s modality MDS stress value (2D): %.4f\n', [upper(modality(1)) modality(2:end)], stress);
        
        % Compute stress for dimensions 1 to 10
        for dim = 1:10
            [~, stress] = mdscale(dissimilarityMatrix, dim, 'Criterion', 'stress');
            stressByDimension(dim, modIdx) = stress;
            fprintf('%s modality MDS stress value (%dD): %.4f\n', [upper(modality(1)) modality(2:end)], dim, stress);
        end
    else
        fprintf('No valid subjects for %s modality.\n', modality);
        avgSimilarityMatrices{modIdx} = [];
        stimuliSets{modIdx} = {};
        mdsConfigurations{modIdx} = [];
        stressByDimension(:, modIdx) = NaN;
    end
end

% Check if both matrices are available
if isempty(avgSimilarityMatrices{1}) || isempty(avgSimilarityMatrices{2})
    error('Cannot compute correlation or distances: missing data for one or both modalities.');
end

% Compute correlation between haptic and visual similarity matrices
hapticMatrix = avgSimilarityMatrices{1};
visualMatrix = avgSimilarityMatrices{2};
n = size(hapticMatrix, 1);
upperTriIdx = triu(true(n), 1); % Upper triangle, excluding diagonal
hapticVector = hapticMatrix(upperTriIdx);
visualVector = visualMatrix(upperTriIdx);
[rho, pval] = corr(hapticVector, visualVector, 'Type', 'Pearson');
fprintf('Pearson correlation between haptic and visual similarity matrices: rho = %.4f, p = %.4f\n', rho, pval);

% Assume two categories based on stimuli labels (e.g., G11-G15 vs. G16-G28)
stimuli = stimuliSets{1}; % Assuming same stimuli for both modalities
category1 = stimuli(contains(stimuli, {'G11', 'G12', 'G13', 'G14', 'G15'}));
category2 = stimuli(~contains(stimuli, {'G11', 'G12', 'G13', 'G14', 'G15'}));
fprintf('Category 1 stimuli: %s\n', strjoin(category1, ', '));
fprintf('Category 2 stimuli: %s\n', strjoin(category2, ', '));

% Compute within- and between-category distances for each modality (no plotting)
withinDistances = zeros(2, 2); % [Haptic, Visual] x [Within1, Within2]
betweenDistances = zeros(1, 2); % [Haptic, Visual]
for modIdx = 1:2
    mdsConfig = mdsConfigurations{modIdx};
    
    % Indices for categories
    cat1Idx = find(ismember(stimuli, category1));
    cat2Idx = find(ismember(stimuli, category2));
    
    % Within-category distances
    distCat1 = pdist(mdsConfig(cat1Idx, :), 'euclidean');
    distCat2 = pdist(mdsConfig(cat2Idx, :), 'euclidean');
    withinDistances(modIdx, 1) = mean(distCat1);
    withinDistances(modIdx, 2) = mean(distCat2);
    
    % Between-category distances
    betweenDist = [];
    for i = 1:length(cat1Idx)
        for j = 1:length(cat2Idx)
            dist = norm(mdsConfig(cat1Idx(i), :) - mdsConfig(cat2Idx(j), :));
            betweenDist = [betweenDist, dist];
        end
    end
    betweenDistances(modIdx) = mean(betweenDist);
    
    fprintf('%s modality: Within-Cat1 = %.4f, Within-Cat2 = %.4f, Between = %.4f\n', ...
        [upper(modality(1)) modality(2:end)], withinDistances(modIdx, 1), withinDistances(modIdx, 2), betweenDistances(modIdx));
end

% Create figure with subplots
figure;

% Haptic: Similarity matrix
subplot(2, 3, 1);
imagesc(avgSimilarityMatrices{1});
colormap(parula);
colorbar;
caxis([1, 7]); % Set consistent color range
title(sprintf('Avg Similarity - Haptic (N=%d)', validSubjectCounts(1)));
xlabel('Stimulus');
ylabel('Stimulus');
set(gca, 'XTick', 1:length(stimuliSets{1}), 'XTickLabel', stimuliSets{1}, 'YTick', 1:length(stimuliSets{1}), 'YTickLabel', stimuliSets{1});
set(gca, 'XTickLabelRotation', 45);
set(gca, 'FontSize', 8);
axis square;

% Haptic: MDS configuration
subplot(2, 3, 2);
scatter(mdsConfigurations{1}(:,1), mdsConfigurations{1}(:,2), 50, 'filled');
hold on;
for i = 1:length(stimuliSets{1})
    text(mdsConfigurations{1}(i,1), mdsConfigurations{1}(i,2), stimuliSets{1}{i}, 'FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
hold off;
title(sprintf('MDS (2D) - Haptic (Stress = %.4f)', stressValues(1)));
xlabel('Dimension 1');
ylabel('Dimension 2');
xlim([-4, 4]); % Set consistent range
ylim([-4, 4]); % Set consistent range
grid on;
set(gca, 'FontSize', 8);
axis square;

% Visual: Similarity matrix
subplot(2, 3, 4);
imagesc(avgSimilarityMatrices{2});
colormap(parula);
colorbar;
caxis([1, 7]); % Set consistent color range
title(sprintf('Avg Similarity - Visual (N=%d)', validSubjectCounts(2)));
xlabel('Stimulus');
ylabel('Stimulus');
set(gca, 'XTick', 1:length(stimuliSets{2}), 'XTickLabel', stimuliSets{2}, 'YTick', 1:length(stimuliSets{2}), 'YTickLabel', stimuliSets{2});
set(gca, 'XTickLabelRotation', 45);
set(gca, 'FontSize', 8);
axis square;

% Visual: MDS configuration
subplot(2, 3, 5);
scatter(mdsConfigurations{2}(:,1), mdsConfigurations{2}(:,2), 50, 'filled');
hold on;
for i = 1:length(stimuliSets{2})
    text(mdsConfigurations{2}(i,1), mdsConfigurations{2}(i,2), stimuliSets{2}{i}, 'FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
hold off;
title(sprintf('MDS (2D) - Visual (Stress = %.4f)', stressValues(2)));
xlabel('Dimension 1');
ylabel('Dimension 2');
xlim([-4, 4]); % Set consistent range
ylim([-4, 4]); % Set consistent range
grid on;
set(gca, 'FontSize', 8);
axis square;

% Stress vs. Dimensions
subplot(2, 3, 3);
plot(1:10, stressByDimension(:,1), '-o', 'LineWidth', 2, 'DisplayName', 'Haptic');
hold on;
plot(1:10, stressByDimension(:,2), '-s', 'LineWidth', 2, 'DisplayName', 'Visual');
hold off;
title('MDS Stress vs. Dimensions');
xlabel('Dimensions');
ylabel('Stress Value');
legend('Location', 'northeast');
grid on;
set(gca, 'FontSize', 8);
set(gca, 'XTick', 1:10);

% Correlation
subplot(2, 3, 6);
scatter(hapticVector, visualVector, 10, 'filled', 'DisplayName', 'Similarity Pairs');
xlabel('Haptic Similarity');
ylabel('Visual Similarity');
title(sprintf('Correlation: rho = %.4f, p = %.4f', rho, pval));
grid on;
set(gca, 'FontSize', 8);
legend('Location', 'northeast');
axis square;

% Adjust figure layout
sgtitle('Similarity, MDS, Stress, and Correlation Analysis');
set(gcf, 'Position', [100, 100, 1800, 800]);
set(gcf, 'Units', 'normalized');
set(gcf, 'OuterPosition', [0 0 1 1]);

% Save the figure
saveas(gcf, 'mds_correlation_analysis.png');