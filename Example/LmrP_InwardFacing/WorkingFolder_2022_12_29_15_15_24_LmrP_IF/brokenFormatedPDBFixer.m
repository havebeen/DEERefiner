function mutatedFormatedPDB = brokenFormatedPDBFixer(formatedPDB)
    % H atom removal 
    HRemovedformatedPDB = HAtomRemover(formatedPDB);

    % find out unknown and broken side chain 
    firstProblematicResiduesIndex = brokenResidueDetector(HRemovedformatedPDB);
    
    % mutated borken and unknown residue to ALA
    firstMutatedFormatedPDB = brokenResidueFixerALA(HRemovedformatedPDB, firstProblematicResiduesIndex);
    
    % find out still unknown and broken side chain
    secondProblematicResiduesIndex = brokenResidueDetector(firstMutatedFormatedPDB);
    
    if any(secondProblematicResiduesIndex) == 0
        mutatedFormatedPDB = firstMutatedFormatedPDB;
    else
        secondMutatedFormatedPDB = brokenResidueFixerGLY(firstMutatedFormatedPDB, secondProblematicResiduesIndex);
        mutatedFormatedPDB = secondMutatedFormatedPDB;
    end
        
end

function HRemovedformatedPDB = HAtomRemover(formatedPDB)
    HRemovedformatedPDB = formatedPDB(formatedPDB(:, 15) ~= " H", :);
end

function problematicResiduesIndex = brokenResidueDetector(formatedPDB)
    alternateRotamerRemovedFormatedPDB = formatedPDB(formatedPDB(:, 4)==" "|formatedPDB(:, 4)=="A", :);
    residueIndex = unique(double(alternateRotamerRemovedFormatedPDB(:, 7)));
    problematicResiduesIndex = zeros(1, length(residueIndex));
    for loopResidueIndex = 1:length(residueIndex)
        residueType = unique(alternateRotamerRemovedFormatedPDB(double(alternateRotamerRemovedFormatedPDB(:, 7))==residueIndex(loopResidueIndex), 5));
        residueTypeMapFileName = strcat(residueType, "_mapping.txt");
        ifMappingFileExist = isfile(residueTypeMapFileName);
        if ifMappingFileExist == 0
            problematicResiduesIndex(loopResidueIndex) = 2;
        end
        residueLength(loopResidueIndex) = length(alternateRotamerRemovedFormatedPDB(double(alternateRotamerRemovedFormatedPDB(:, 7))==residueIndex(loopResidueIndex), 1));
        if problematicResiduesIndex(loopResidueIndex) == 0
            residueTypeMapFileIndex = fopen(residueTypeMapFileName, 'r');
            residueTypeMapCell = textscan(residueTypeMapFileIndex, '%s\t', 'CommentStyle', '#');
            fclose(residueTypeMapFileIndex);
            residueTypeMapStringMat = string(residueTypeMapCell{1});   % transform cell to string
            residueTypeMapStringMat(1) = []; % remove residue name
            residueTypeMapStringMat = reshape(residueTypeMapStringMat, 2, [])';   % reshape string matrix
            residueTypeSideChainAtom = residueTypeMapStringMat(:, 1);
            if length(residueTypeSideChainAtom) ~= residueLength(loopResidueIndex)
                problematicResiduesIndex(loopResidueIndex) = 1;
            end
        end
    end
end

function mutatedFormatedPDB = brokenResidueFixerALA(formatedPDB, problematicResiduesIndex)
    alternateRotamerRemovedFormatedPDB = formatedPDB(formatedPDB(:, 4)==" "|formatedPDB(:, 4)=="A", :);
    residueIndex = unique(double(alternateRotamerRemovedFormatedPDB(:, 7)));
    for loopResidueIndex = 1:length(residueIndex)
        if problematicResiduesIndex(loopResidueIndex) ~=0
            alternateRotamerRemovedFormatedPDB(double(alternateRotamerRemovedFormatedPDB(:, 7))==residueIndex(loopResidueIndex), 5) = "ALA";
            alternateRotamerRemovedFormatedPDB(double(alternateRotamerRemovedFormatedPDB(:, 7))==residueIndex(loopResidueIndex)&...
            alternateRotamerRemovedFormatedPDB(:, 3)~=" N  "&...
            alternateRotamerRemovedFormatedPDB(:, 3)~=" CA "&...
            alternateRotamerRemovedFormatedPDB(:, 3)~=" C  "&...
            alternateRotamerRemovedFormatedPDB(:, 3)~=" O  "&...
            alternateRotamerRemovedFormatedPDB(:, 3)~=" CB ", :) = [];
        end
       
    end
    mutatedFormatedPDB = alternateRotamerRemovedFormatedPDB;
end

function mutatedFormatedPDB = brokenResidueFixerGLY(formatedPDB, problematicResiduesIndex)
    alternateRotamerRemovedFormatedPDB = formatedPDB(formatedPDB(:, 4)==" "|formatedPDB(:, 4)=="A", :);
    residueIndex = unique(double(alternateRotamerRemovedFormatedPDB(:, 7)));
    for loopResidueIndex = 1:length(residueIndex)
        if problematicResiduesIndex(loopResidueIndex) ~=0
            alternateRotamerRemovedFormatedPDB(double(alternateRotamerRemovedFormatedPDB(:, 7))==residueIndex(loopResidueIndex), 5) = "GLY";
            alternateRotamerRemovedFormatedPDB(double(alternateRotamerRemovedFormatedPDB(:, 7))==residueIndex(loopResidueIndex)&...
            alternateRotamerRemovedFormatedPDB(:, 3)~=" N  "&...
            alternateRotamerRemovedFormatedPDB(:, 3)~=" CA "&...
            alternateRotamerRemovedFormatedPDB(:, 3)~=" C  "&...
            alternateRotamerRemovedFormatedPDB(:, 3)~=" O  ", :) = [];
        end
       
    end
    mutatedFormatedPDB = alternateRotamerRemovedFormatedPDB;
end