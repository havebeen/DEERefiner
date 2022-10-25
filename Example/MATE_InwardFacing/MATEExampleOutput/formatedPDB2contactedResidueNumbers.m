function contactedResidueIndex = formatedPDB2contactedResidueNumbers(formatedPDB, contactCriterion)
    formatedPDBCoordinates = double(formatedPDB(:, 9:11));
    pairDistances = pdist2(formatedPDBCoordinates, formatedPDBCoordinates);
    residueIndex = unique(double(formatedPDB(:, 7)));
    
    selfIgnoringMatrix = diag(NaN(length(formatedPDBCoordinates(:, 1)), 1));
    ignoredPairDistance = pairDistances+selfIgnoringMatrix;
    
    for i = 1:length(residueIndex)
        ignoredPairDistance(double(formatedPDB(:, 7)) == (residueIndex(i)) & (formatedPDB(:, 3) == " N  "|formatedPDB(:, 3) == " CA "|formatedPDB(:, 3) == " C  "|formatedPDB(:, 3) == " O  "|formatedPDB(:, 3) == " CB "),...
            double(formatedPDB(:, 7)) == (residueIndex(i)) & (formatedPDB(:, 3) == " N  "|formatedPDB(:, 3) == " CA "|formatedPDB(:, 3) == " C  "|formatedPDB(:, 3) == " O  "|formatedPDB(:, 3) == " CB ")) = NaN;
        ignoredPairDistance(double(formatedPDB(:, 7)) == (residueIndex(i)) &((formatedPDB(:, 3) ~= " N  "&formatedPDB(:, 3) ~= " CA "&formatedPDB(:, 3) ~= " C  "&formatedPDB(:, 3) ~= " O  "&formatedPDB(:, 3) ~= " CB ")),...
            double(formatedPDB(:, 7)) == (residueIndex(i)) &((formatedPDB(:, 3) ~= " N  "&formatedPDB(:, 3) ~= " CA "&formatedPDB(:, 3) ~= " C  "&formatedPDB(:, 3) ~= " O  "&formatedPDB(:, 3) ~= " CB "))) = NaN;
        ignoredPairDistance(double(formatedPDB(:, 7)) == (residueIndex(i)) & formatedPDB(:, 3) == " CB ",...
            double(formatedPDB(:, 7)) == (residueIndex(i)) & (formatedPDB(:, 3) == " CG1" |formatedPDB(:, 3) == " CG2"|formatedPDB(:, 3) == " CG "|formatedPDB(:, 3) == " OG "|formatedPDB(:, 3) == " OG1"|formatedPDB(:, 3) == " SG ")) = NaN;
        ignoredPairDistance(double(formatedPDB(:, 7)) == (residueIndex(i)) & (formatedPDB(:, 5) == "PRO"), ...
            double(formatedPDB(:, 7)) == (residueIndex(i)) & (formatedPDB(:, 5) == "PRO")) = NaN;
        ignoredPairDistance(double(formatedPDB(:, 7)) == (residueIndex(i)) & (formatedPDB(:, 5) == "ALA"), ...
            double(formatedPDB(:, 7)) == (residueIndex(i)) & (formatedPDB(:, 5) == "ALA")) = NaN;
        ignoredPairDistance(double(formatedPDB(:, 7)) == (residueIndex(i)) & (formatedPDB(:, 5) == "GLY"), ...
            double(formatedPDB(:, 7)) == (residueIndex(i)) & (formatedPDB(:, 5) == "GLY")) = NaN;
    end
    
    for i = 1:length(residueIndex)-1
        ignoredPairDistance(double(formatedPDB(:, 7)) == (residueIndex(i)) & formatedPDB(:, 3) == " C  ", double(formatedPDB(:, 7)) == (residueIndex(i+1)) & (formatedPDB(:, 3) == " N  ")) = NaN;
        ignoredPairDistance(double(formatedPDB(:, 7)) == (residueIndex(i)) & formatedPDB(:, 3) == " O  ", double(formatedPDB(:, 7)) == (residueIndex(i+1)) & (formatedPDB(:, 3) == " N  ")) = NaN;
    end
    
    ignoredFinalPairDistance = (ignoredPairDistance+ignoredPairDistance')/2;
    clashedSideChain = formatedPDB(sum(ignoredFinalPairDistance<contactCriterion)>0, :);
    clashedSideChain(clashedSideChain(:, 3) == " N  "|clashedSideChain(:, 3) == " CA "|clashedSideChain(:, 3) == " C  "|clashedSideChain(:, 3) == " O  "|clashedSideChain(:, 3) == " CB ", :) = [];
    contactedResidueIndex = unique(double(clashedSideChain(:, 7)));
    
end