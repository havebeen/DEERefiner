function [residue1List, residue2List, distanceDistributionList] = advancedPrFileLoader(advancedPrFileName)
    advancedPrFileIndex = fopen(advancedPrFileName, 'r');
    
    DEERDistanceDistributionFileNameCell= textscan(advancedPrFileIndex, '%s\n');
    
    fclose(advancedPrFileIndex);
    
    DEERDistanceDistributionFileName = string(DEERDistanceDistributionFileNameCell{1});
    
    for DEERDistanceDistributionFileIndex = 1:length(DEERDistanceDistributionFileName)
        
        [residue1List(DEERDistanceDistributionFileIndex), residue2List(DEERDistanceDistributionFileIndex), distanceDistributionList(:, :, DEERDistanceDistributionFileIndex)] = PrFileReader(DEERDistanceDistributionFileName(DEERDistanceDistributionFileIndex));
        
        distanceDistributionList(:, 2, DEERDistanceDistributionFileIndex) = distanceDistributionList(:, 2, DEERDistanceDistributionFileIndex)/sum(distanceDistributionList(:, 2, DEERDistanceDistributionFileIndex));
        
    end
    
end