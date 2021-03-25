for n = 1: ms.numNeurons;
FiringFields(n,:) = extended_calcium_analysis_output(n).likelihood';
end

[~,maxFiringFields] = max(FiringFields, [], 2);
[~, index] = sort(maxFiringFields, 'ascend');
sortedFiringFields = FiringFields(index, :);

clims = [0 0.15];
imagesc(sortedFiringFields, clims);
colorbar