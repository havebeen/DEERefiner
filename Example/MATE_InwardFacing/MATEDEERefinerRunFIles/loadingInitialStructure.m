function [initialStructureFullPath, initialStructurePathName, initialStructureFileName] = loadingInitialStructure
    [initialStructureFileName, initialStructurePathName] = uigetfile('*.pdb');
    if initialStructureFileName==0
        initialStructureFullPath = "";
      return
    end
    cd(initialStructurePathName)
    initialStructureFullPath = strcat(initialStructurePathName, initialStructureFileName);
end