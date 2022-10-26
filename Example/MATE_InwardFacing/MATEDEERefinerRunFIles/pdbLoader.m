function formatedPDB = pdbLoader(targetPDBFileName)
    
    % generating PDB string cell
    targetPDBFileStringCell = targetPDBFileStringCellGenerator(targetPDBFileName);
    
    % generating PDB char array
    targetPDBFileChar = targetPDBFileCharGenerator(targetPDBFileStringCell);
    
    % fixxing MD PDB file
    targetPDBFileChar = MDPDBFixer(targetPDBFileChar);
    
    % checking if the input is a NMR structure
    NMRCheckedFormatedPDBchar = NMRFinalCharGenerator(targetPDBFileChar);
    
    % chekcing if the input is multi-chains crystal structure
    finalFormatedPDBchar = XrayFinalCharGenerator(NMRCheckedFormatedPDBchar);
    
    % generating formated PDB string array
    formatedPDB = atomPartFormator(finalFormatedPDBchar);

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

function NMRCheckedFormatedPDBchar = NMRFinalCharGenerator(targetPDBFileChar)
    targetPDBMODELPart = targetPDBFileChar(all(targetPDBFileChar(:, 1:5)=='MODEL', 2), :);
    targetPDBMODELLineNumber = length(targetPDBMODELPart(:, 1));
    if targetPDBMODELLineNumber>1
        NMRStructureModelLineIndex = find(all(ismember(targetPDBFileChar(:, 1:5), 'MODEL'), 2));
        NMRCheckedFormatedPDBchar = targetPDBFileChar(NMRStructureModelLineIndex(1)+1:NMRStructureModelLineIndex(2)-1, :);
    else
        NMRCheckedFormatedPDBchar = targetPDBFileChar;
    end
end

function finalFormatedPDBchar = XrayFinalCharGenerator(NMRCheckedFormatedPDBchar)
    targetPDBATOMPart = NMRCheckedFormatedPDBchar(all(NMRCheckedFormatedPDBchar(:, 1:4)=='ATOM', 2), :);
    finalFormatedPDBchar = targetPDBATOMPart(targetPDBATOMPart(:, 22) == 'A', :);
end

function formatedAtomPart = atomPartFormator(targetPDBFileChar)
    targetPDBATOMPart = targetPDBFileChar(all(targetPDBFileChar(:, 1:4)=='ATOM', 2), :);
    formatedAtomPart(:, 1) = string(targetPDBATOMPart(:, 1:4));
    formatedAtomPart(:, 2) = string(targetPDBATOMPart(:, 7:11));
    formatedAtomPart(:, 3) = string(targetPDBATOMPart(:, 13:16));
    formatedAtomPart(:, 4) = string(targetPDBATOMPart(:, 17));
    formatedAtomPart(:, 5) = string(targetPDBATOMPart(:, 18:20));
    formatedAtomPart(:, 6) = string(targetPDBATOMPart(:, 22));
    formatedAtomPart(:, 7) = string(targetPDBATOMPart(:, 23:26));
    formatedAtomPart(:, 8) = string(targetPDBATOMPart(:, 27));
    formatedAtomPart(:, 9) = string(targetPDBATOMPart(:, 31:38));
    formatedAtomPart(:, 10) = string(targetPDBATOMPart(:, 39:46));
    formatedAtomPart(:, 11) = string(targetPDBATOMPart(:, 47:54));
    formatedAtomPart(:, 12) = string(targetPDBATOMPart(:, 55:60));
    formatedAtomPart(:, 13) = string(targetPDBATOMPart(:, 61:66));
    formatedAtomPart(:, 14) = string(targetPDBATOMPart(:, 73:76));
    formatedAtomPart(:, 15) = string(targetPDBATOMPart(:, 77:78));
    formatedAtomPart(:, 16) = string(targetPDBATOMPart(:, 79:80));
end

















































