function contactedResidueIndex = formatedPDB2contactedResidueNumbers(formatedPDB, contactCriterion)
    formatedPDB(:, 15) = "  1.00";
    formatedPDB(formatedPDB(:, 3) == " N  "|formatedPDB(:, 3) == " CA "|formatedPDB(:, 3) == " C  "|formatedPDB(:, 3) == " O  "|formatedPDB(:, 3) == " CB ", 15) = "  0.00";
    
    formatedPDBCoordinates = double(formatedPDB(:, 9:11));
    pairDistances = pdist2(formatedPDBCoordinates, formatedPDBCoordinates);
    residueIndex = unique(double(formatedPDB(:, 7)));
    
    selfIgnoringMatrix = diag(NaN(length(formatedPDBCoordinates(:, 1)), 1));
    ignoredPairDistance = pairDistances+selfIgnoringMatrix;
    
    residueIndexString = strrep(formatedPDB(:, 7), ' ', '');
    for i = 1:length(residueIndex)
        ignoredPairDistance(residueIndexString == string(residueIndex(i)) & (formatedPDB(:, 15) == "  0.00"),...
            residueIndexString == string(residueIndex(i)) & (formatedPDB(:, 15) == "  0.00")) = NaN;
        ignoredPairDistance(residueIndexString == string(residueIndex(i)) & (formatedPDB(:, 15) == "  1.00"),...
            residueIndexString == string(residueIndex(i)) & (formatedPDB(:, 15) == "  1.00")) = NaN;
        ignoredPairDistance(residueIndexString == string(residueIndex(i)) & formatedPDB(:, 3) == " CA ",...
            residueIndexString == string(residueIndex(i)) & contains(formatedPDB(:, 3), 'G')) = NaN;
        ignoredPairDistance(residueIndexString == string(residueIndex(i)) & formatedPDB(:, 3) == " CB ",...
            residueIndexString == string(residueIndex(i)) & contains(formatedPDB(:, 3), 'D')) = NaN;
        ignoredPairDistance(residueIndexString == string(residueIndex(i)) & formatedPDB(:, 3) == " CB ",...
            residueIndexString == string(residueIndex(i)) & contains(formatedPDB(:, 3), 'G')) = NaN;
        
%         ignoredPairDistance(residueIndexString == string(residueIndex(i)) & (formatedPDB(:, 3) == " N  "|formatedPDB(:, 3) == " CA "|formatedPDB(:, 3) == " C  "|formatedPDB(:, 3) == " O  "|formatedPDB(:, 3) == " CB "),...
%             residueIndexString == string(residueIndex(i)) & (formatedPDB(:, 3) == " N  "|formatedPDB(:, 3) == " CA "|formatedPDB(:, 3) == " C  "|formatedPDB(:, 3) == " O  "|formatedPDB(:, 3) == " CB ")) = NaN;
%         ignoredPairDistance(residueIndexString == string(residueIndex(i)) &((formatedPDB(:, 3) ~= " N  "&formatedPDB(:, 3) ~= " CA "&formatedPDB(:, 3) ~= " C  "&formatedPDB(:, 3) ~= " O  "&formatedPDB(:, 3) ~= " CB ")),...
%             residueIndexString == string(residueIndex(i)) &((formatedPDB(:, 3) ~= " N  "&formatedPDB(:, 3) ~= " CA "&formatedPDB(:, 3) ~= " C  "&formatedPDB(:, 3) ~= " O  "&formatedPDB(:, 3) ~= " CB "))) = NaN;
%         ignoredPairDistance(residueIndexString == string(residueIndex(i)) & formatedPDB(:, 3) == " CB ",...
%             residueIndexString == string(residueIndex(i)) & (formatedPDB(:, 3) == " CG1" |formatedPDB(:, 3) == " CG2"|formatedPDB(:, 3) == " CG "|formatedPDB(:, 3) == " OG "|formatedPDB(:, 3) == " OG1"|formatedPDB(:, 3) == " SG "|...
%                                                               formatedPDB(:, 3) == " CD " |formatedPDB(:, 3) == " CD1"|formatedPDB(:, 3) == " CD2"|formatedPDB(:, 3) == " OD1"|formatedPDB(:, 3) == " OD2"|formatedPDB(:, 3) == " SD "|formatedPDB(:, 3) == " ND1"|formatedPDB(:, 3) == " ND2")) = NaN;
%         ignoredPairDistance(residueIndexString == string(residueIndex(i)) & formatedPDB(:, 3) == " CA ",...
%             residueIndexString == string(residueIndex(i)) & (formatedPDB(:, 3) == " CG1" |formatedPDB(:, 3) == " CG2"|formatedPDB(:, 3) == " CG "|formatedPDB(:, 3) == " OG "|formatedPDB(:, 3) == " OG1"|formatedPDB(:, 3) == " SG ")) = NaN;
        ignoredPairDistance(residueIndexString == string(residueIndex(i)) & (formatedPDB(:, 5) == "PRO"), ...
            residueIndexString == string(residueIndex(i)) & (formatedPDB(:, 5) == "PRO")) = NaN;
        ignoredPairDistance(residueIndexString == string(residueIndex(i)) & (formatedPDB(:, 5) == "ALA"), ...
            residueIndexString == string(residueIndex(i)) & (formatedPDB(:, 5) == "ALA")) = NaN;
        ignoredPairDistance(residueIndexString == string(residueIndex(i)) & (formatedPDB(:, 5) == "GLY"), ...
            residueIndexString == string(residueIndex(i)) & (formatedPDB(:, 5) == "GLY")) = NaN;
        ignoredPairDistance(residueIndexString == string(residueIndex(i)) & formatedPDB(:, 3) == " O  ", :) = NaN;
    end
    
    for i = 1:length(residueIndex)-1
        ignoredPairDistance(residueIndexString == string(residueIndex(i)) & formatedPDB(:, 3) == " C  ", residueIndexString == string(residueIndex(i+1)) & (formatedPDB(:, 3) == " N  ")) = NaN;
        ignoredPairDistance(residueIndexString == string(residueIndex(i)) & formatedPDB(:, 3) == " O  ", residueIndexString == string(residueIndex(i+1)) & (formatedPDB(:, 3) == " N  ")) = NaN;
    end
    
    ignoredFinalPairDistance = (ignoredPairDistance+ignoredPairDistance')/2;
    clashedSideChain = formatedPDB(sum(ignoredFinalPairDistance<contactCriterion)>0, :);
    clashedSideChain(clashedSideChain(:, 3) == " N  "|clashedSideChain(:, 3) == " CA "|clashedSideChain(:, 3) == " C  "|clashedSideChain(:, 3) == " O  "|clashedSideChain(:, 3) == " CB ", :) = [];
    contactedResidueIndex = unique(double(clashedSideChain(:, 7)));
    
end