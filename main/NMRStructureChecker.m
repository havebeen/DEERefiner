function NMRStructureChecker(app)
    inputStructure = app.initialStructureFileName;
    targetPDBFileStringCell = targetPDBFileStringCellGenerator(inputStructure);
    targetPDBFileChar = targetPDBFileCharGenerator(targetPDBFileStringCell);
    targetPDBFileChar = MDPDBFixer(targetPDBFileChar);
    targetPDBMODELPart = targetPDBFileChar(all(targetPDBFileChar(:, 1:5)=='MODEL', 2), :);
    targetPDBMODELLineNumber = length(targetPDBMODELPart(:, 1));
    if targetPDBMODELLineNumber>1
        
        NMREnsembleFormatedPDBchar = NMREnsembleCharGenerator(targetPDBFileChar);
        formatedAtomPart = atomPartFormator(targetPDBFileChar, NMREnsembleFormatedPDBchar);


        timeStamp = char(erase(app.runFileName(1),"runFile_"));
        timeStamp(end-1:end) = [];
        advancedPrFileName = strcat(timeStamp, '_advanced_Pr.dat');
        [residue1List, residue2List, distanceDistributionList] = advancedPrFileLoader(advancedPrFileName);
        [~, targetDEERMaximalIndexTmp] = max(distanceDistributionList);
    
        targetDEERMaximalIndex = reshape(targetDEERMaximalIndexTmp(1, 2, :), 1, []);
    
        R1LibraryMatFileName = 'R1_library_20210523.mat';
        load(R1LibraryMatFileName);
        [~, targetDistanceDistributionFinal] = targetDistanceDistributionFinalForJSDGenerator(app.distanceDistributionFullPath);
        [distanceDistributionFinal, ~] = distanceDistributionFinalGenerator_internal(formatedAtomPart, residue1List, residue2List, R1_20210523);
        for FINALPDBNumber = 1:length(formatedAtomPart(1, 1, :))
            [~, FINALPDBMaximalIndex(FINALPDBNumber, :)] = max(distanceDistributionFinal(:, :, FINALPDBNumber));
        end
        differenceBetweenFINALPDBAndTarget = abs(FINALPDBMaximalIndex-targetDEERMaximalIndex);
        FINALPDBMaximalDifference = max(differenceBetweenFINALPDBAndTarget');
        FINALPDBRMSE = sqrt(mean((differenceBetweenFINALPDBAndTarget.^2)'));
        minimalFINALPDBMaximalDifference = min(FINALPDBMaximalDifference);
        minimalFINALPDBRMSE = min(FINALPDBRMSE(FINALPDBMaximalDifference == min(FINALPDBMaximalDifference)));
        minimalRMSEIndex = minimalRMSEPDBGenerator_internal(FINALPDBMaximalDifference, minimalFINALPDBMaximalDifference, FINALPDBRMSE, minimalFINALPDBRMSE, timeStamp);
        newFileName = strcat(erase(app.initialStructureFileName, '.pdb'), '_', num2str(minimalRMSEIndex), '.pdb');
        pdbSaver(newFileName, formatedAtomPart(:, :, minimalRMSEIndex))
        app.initialStructureFileName = newFileName;
        app.initialStructureFullPath = strcat(app.initialStructurePath, app.initialStructureFileName);
    else
        return
    end
end





function targetPDBFileStringCell = targetPDBFileStringCellGenerator(targetPDBFileName)
    targetPDBFileIndex = fopen(targetPDBFileName, 'r');
    targetPDBFileStringCell = textscan(targetPDBFileIndex,'%s', 'delimiter', '\n');
    fclose(targetPDBFileIndex);
end


function targetPDBFileChar = targetPDBFileCharGenerator(targetPDBFileStringCell)
    targetPDBFileChar = char(targetPDBFileStringCell{1});
end

function targetPDBFileChar = MDPDBFixer(targetPDBFileChar)
    targetPDBFileLineLength = length(targetPDBFileChar(1, :));
    if targetPDBFileLineLength == 78
        targetPDBFileChar = [targetPDBFileChar repmat('  ', length(targetPDBFileChar(:, 1)), 1)];
    end
end

function NMREnsembleFormatedPDBchar = NMREnsembleCharGenerator(targetPDBFileChar)
    targetPDBMODELPart = targetPDBFileChar(all(targetPDBFileChar(:, 1:5)=='MODEL', 2), :);
    targetPDBMODELLineNumber = length(targetPDBMODELPart(:, 1));
    
    
    if targetPDBMODELLineNumber>1
        targetPDBATOMPart = targetPDBFileChar(all(targetPDBFileChar(:, 1:4)=='ATOM', 2), :);
        singleStructureLength = length(targetPDBATOMPart)/targetPDBMODELLineNumber;
        NMREnsembleFormatedPDBchar = pagetranspose(reshape(targetPDBATOMPart', 80, singleStructureLength, []));
    else
        NMREnsembleFormatedPDBchar = targetPDBFileChar;
    end
end

function formatedAtomPart = atomPartFormator(targetPDBFileChar, NMREnsembleFormatedPDBchar)
    targetPDBMODELPart = targetPDBFileChar(all(targetPDBFileChar(:, 1:5)=='MODEL', 2), :);
    targetPDBMODELLineNumber = length(targetPDBMODELPart(:, 1));
    for NMRStructureIndex = 1:targetPDBMODELLineNumber
    formatedAtomPart(:, 1, NMRStructureIndex) = string(NMREnsembleFormatedPDBchar(:, 1:4, NMRStructureIndex));
    formatedAtomPart(:, 2, NMRStructureIndex) = string(NMREnsembleFormatedPDBchar(:, 7:11, NMRStructureIndex));
    formatedAtomPart(:, 3, NMRStructureIndex) = string(NMREnsembleFormatedPDBchar(:, 13:16, NMRStructureIndex));
    formatedAtomPart(:, 4, NMRStructureIndex) = string(NMREnsembleFormatedPDBchar(:, 17, NMRStructureIndex));
    formatedAtomPart(:, 5, NMRStructureIndex) = string(NMREnsembleFormatedPDBchar(:, 18:20, NMRStructureIndex));
    formatedAtomPart(:, 6, NMRStructureIndex) = string(NMREnsembleFormatedPDBchar(:, 22, NMRStructureIndex));
    formatedAtomPart(:, 7, NMRStructureIndex) = string(NMREnsembleFormatedPDBchar(:, 23:26, NMRStructureIndex));
    formatedAtomPart(:, 8, NMRStructureIndex) = string(NMREnsembleFormatedPDBchar(:, 27, NMRStructureIndex));
    formatedAtomPart(:, 9, NMRStructureIndex) = string(NMREnsembleFormatedPDBchar(:, 31:38, NMRStructureIndex));
    formatedAtomPart(:, 10, NMRStructureIndex) = string(NMREnsembleFormatedPDBchar(:, 39:46, NMRStructureIndex));
    formatedAtomPart(:, 11, NMRStructureIndex) = string(NMREnsembleFormatedPDBchar(:, 47:54, NMRStructureIndex));
    formatedAtomPart(:, 12, NMRStructureIndex) = string(NMREnsembleFormatedPDBchar(:, 55:60, NMRStructureIndex));
    formatedAtomPart(:, 13, NMRStructureIndex) = string(NMREnsembleFormatedPDBchar(:, 61:66, NMRStructureIndex));
    formatedAtomPart(:, 14, NMRStructureIndex) = string(NMREnsembleFormatedPDBchar(:, 73:76, NMRStructureIndex));
    formatedAtomPart(:, 15, NMRStructureIndex) = string(NMREnsembleFormatedPDBchar(:, 77:78, NMRStructureIndex));
    formatedAtomPart(:, 16, NMRStructureIndex) = string(NMREnsembleFormatedPDBchar(:, 79:80, NMRStructureIndex));
    end
end

function [distanceDistributionFinal, FINALPDBENSEMBLE] = distanceDistributionFinalGenerator_internal(FINALPDBENSEMBLE, residue1List, residue2List, R1_20210523)
%     FINALPDBENSEMBLE = pdbModelsLoader(FINALPDBENSEMBLEFileName);
    for FINALPDBNumber = 1:length(FINALPDBENSEMBLE(1, 1, :))
        for labelingIndex = 1:length(residue1List)
            [distanceDistributionFinal(:, labelingIndex, FINALPDBNumber), ~] = DEERefineMTSSLLabeling(FINALPDBENSEMBLE(:, :, FINALPDBNumber), residue1List(labelingIndex), residue2List(labelingIndex), R1_20210523);
        end
    end
end


function minimalRMSEIndex = minimalRMSEPDBGenerator_internal(FINALPDBMaximalDifference, minimalFINALPDBMaximalDifference, FINALPDBRMSE, minimalFINALPDBRMSE, timeStamp)
    FINALPDBFileName = strcat("FINALPDB_",  timeStamp, ".pdb");
    minimalRMSEIndex = find((FINALPDBMaximalDifference == minimalFINALPDBMaximalDifference &FINALPDBRMSE== minimalFINALPDBRMSE)==1);
%     pdbSaver(FINALPDBFileName, FINALPDBENSEMBLE(:, :, minimalRMSEIndex(1)));
end



