function DEERefineFileCopier(app)
    
    copyfile(app.initialStructureFullPath, '.');
    for distanceDistributionNumber = 1:length(app.distanceDistributionFullPath)
        copyfile(app.distanceDistributionFullPath(distanceDistributionNumber), '.');
    end
    
    load('DEERefineFilesName.mat');
    for DEERefineFilesNumber = 1:length(DEERefineFilesName)
        copyFileFullPathTmp = which(DEERefineFilesName(DEERefineFilesNumber));
        copyfile(copyFileFullPathTmp, '.');
        clear copyFileFullPathTmp
    end
    
end