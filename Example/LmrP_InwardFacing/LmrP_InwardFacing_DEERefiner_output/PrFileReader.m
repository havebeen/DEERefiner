function [residue1, residue2, distanceDistribution] = PrFileReader(PrFileName)
PrFileIndex = fopen(PrFileName, 'r');
distanceDistributionCell = textscan(PrFileIndex, '%f\t%f\n', 'CommentStyle', '#');
fclose(PrFileIndex);

distanceDistributionXAxis = distanceDistributionCell{1};
distanceDistributionYAxis = distanceDistributionCell{2}/max(distanceDistributionCell{2});
distanceDistribution = [distanceDistributionXAxis distanceDistributionYAxis];

PrFileIndex = fopen(PrFileName, 'r');
residueCell = textscan(PrFileIndex, '%s\n', 2);
fclose(PrFileIndex);
residue1 = double(erase(string(residueCell{1}(1)), '#'));
residue2 = double(erase(string(residueCell{1}(2)), '#'));       