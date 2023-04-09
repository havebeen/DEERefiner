function [simulatedDistanceDistributionY, simulatedDistanceDistributionX] = DEERefineMTSSLLabeling(formatedPDB, residue1, residue2, R1_20210523)

    [~, formatedR1AtResidue1] = R1_searching(formatedPDB, R1_20210523, residue1);
    [~, formatedR1AtResidue2] = R1_searching(formatedPDB, R1_20210523, residue2);
    [simulatedDistanceDistributionY, simulatedDistanceDistributionX] = histogramCalculator(formatedR1AtResidue1, formatedR1AtResidue2);
    
end



function [availableR1Index, formatedR1] = R1_searching(formatedPDB, R1_20210523, residueIndex)
    residueIndexString = strrep(formatedPDB(:, 7), ' ', '');
    proteinCoordinates = double(formatedPDB(:, 9:11));
    targetCAlphaCoordinates = double(formatedPDB(formatedPDB(:, 3) == " CA " & residueIndexString == string(residueIndex), 9:11));
    proteinCoordinatesFirstMoved = proteinCoordinates-targetCAlphaCoordinates;
    targetCCoordinatesFirstMoved = proteinCoordinatesFirstMoved(formatedPDB(:, 2) == (formatedPDB(formatedPDB(:, 3) == " C  " & residueIndexString == string(residueIndex), 2)), :);
    vectorCAlphaCFirstMoved = targetCCoordinatesFirstMoved;
    crossCAlphaCZFirstMoved = cross(vectorCAlphaCFirstMoved, [0 0 1]);
    angleCAlphaCZFirstMoved = acos(dot(vectorCAlphaCFirstMoved, [0 0 1])/norm(vectorCAlphaCFirstMoved)/norm([0 0 1]));
    proteinCoordinatesSecondMoved = proteinCoordinatesFirstMoved*rotatingMatrixGenerator(crossCAlphaCZFirstMoved, angleCAlphaCZFirstMoved);
    targetNCoordinatesSecondMoved = proteinCoordinatesSecondMoved(formatedPDB(:, 2) == (formatedPDB(formatedPDB(:, 3) == " N  " & residueIndexString == string(residueIndex), 2)), :);
    vectorCAlphaNSecondMoved = targetNCoordinatesSecondMoved;
    projectedVectorCAlphaNSecondMoved = [vectorCAlphaNSecondMoved(1:2) 0];
    rotatingSignSecondMoved = -1*sign(dot(cross(projectedVectorCAlphaNSecondMoved, [1 0 0]), proteinCoordinatesSecondMoved(formatedPDB(:, 2) == (formatedPDB(formatedPDB(:, 3) == " N  " & residueIndexString == string(residueIndex), 2)), :)));
    angleCAlphaNXSecondMoved = acos(dot(projectedVectorCAlphaNSecondMoved, [1 0 0])/norm(projectedVectorCAlphaNSecondMoved)/norm([1 0 0]));
    proteinCoordinatesThirdMoved = proteinCoordinatesSecondMoved*rotatingMatrixGenerator(proteinCoordinatesSecondMoved(formatedPDB(:, 2) == (formatedPDB(formatedPDB(:, 3) == " C  " & residueIndexString == string(residueIndex), 2)), :), rotatingSignSecondMoved.*angleCAlphaNXSecondMoved);
    rotatedFormatedPDB = formatedPDB;
    rotatedFormatedPDB(:, 9:11) = string(proteinCoordinatesThirdMoved);
    rotatedResidueIndexString = strrep(rotatedFormatedPDB(:, 7), ' ', '');
    
    R1LibraryCoordinates = reshape(pagetranspose(double(R1_20210523(:, 7:9, :))), 3, [])';
    maxMinX = [max(R1LibraryCoordinates(:,1))+1 min(R1LibraryCoordinates(:,1))-1];
    maxMinY = [max(R1LibraryCoordinates(:,2))+1 min(R1LibraryCoordinates(:,2))-1];
    maxMinZ = [max(R1LibraryCoordinates(:,3))+1 min(R1LibraryCoordinates(:,3))-1];
    newProteinCoordinates = double(rotatedFormatedPDB(:, 9:11));
    allAtomIndex = double(rotatedFormatedPDB(rotatedResidueIndexString == string(residueIndex), 2));
    clashedCandidatesAtomIndex = setdiff(double(rotatedFormatedPDB((newProteinCoordinates(:, 1) <= maxMinX(1) &...
                                                                  newProteinCoordinates(:, 1) >= maxMinX(2) &...
                                                                  newProteinCoordinates(:, 2) <= maxMinY(1) &...
                                                                  newProteinCoordinates(:, 2) >= maxMinY(2) &...
                                                                  newProteinCoordinates(:, 3) <= maxMinZ(1) &...
                                                                  newProteinCoordinates(:, 3) >= maxMinZ(2)), 2)), ...
                                                                  allAtomIndex);
    clashedCandidatesAtomCoordinates = newProteinCoordinates(ismember(double(rotatedFormatedPDB(:, 2)),clashedCandidatesAtomIndex), :);
    pairDistanceBetweenClashedCandidatesAndR1Library = pdist2(clashedCandidatesAtomCoordinates, R1LibraryCoordinates);
    R1Number = sum(R1LibraryCoordinates(:,1) == 0 &...
                R1LibraryCoordinates(:,2) == 0 &...
                R1LibraryCoordinates(:,3) == 0);
    R1Index = 1:R1Number;   
    R1Length = length(R1LibraryCoordinates(:,1))/R1Number;
    
    reshapedPairDistanceBetweenClashedCandidatesAndR1Library = reshape(pairDistanceBetweenClashedCandidatesAndR1Library, length(pairDistanceBetweenClashedCandidatesAndR1Library(:, 1)), R1Length, R1Number);
    clashedR1Candidates = reshapedPairDistanceBetweenClashedCandidatesAndR1Library(:, 7:R1Length, :);
    minimalR1Distance = reshape(min(min(clashedR1Candidates)), R1Number, 1);
    [lessConstraintR1Value200, lessConstraintR1Index200] = maxk(minimalR1Distance, 200);
    availableR1Index = lessConstraintR1Index200(lessConstraintR1Value200>1.7);
    availableR1Coordinates = double(R1_20210523(:, 7:9, availableR1Index));
    reshapedAvailableR1Coordinates = reshape(pagetranspose(availableR1Coordinates), 3, [])';
    
    R1FirstMoved = reshapedAvailableR1Coordinates*rotatingMatrixGenerator(proteinCoordinatesThirdMoved(formatedPDB(:, 2) == (formatedPDB(formatedPDB(:, 3) == " C  " & residueIndexString == string(residueIndex), 2)), :), -1*rotatingSignSecondMoved.*angleCAlphaNXSecondMoved);
    R1SecondMoved = R1FirstMoved*rotatingMatrixGenerator(crossCAlphaCZFirstMoved, -1*angleCAlphaCZFirstMoved);
    R1ThirdMoved = R1SecondMoved + targetCAlphaCoordinates;
    reshapedR1ThirdMoved = pagectranspose(reshape(R1ThirdMoved', 3, 18, []));
    formatedR1 = R1_20210523(:, :, availableR1Index);
    formatedR1(:, 7:9, :) = string(reshapedR1ThirdMoved);
end



function [normalizedSimulatedDistanceDistributionY, simulatedDistanceDistributionX] = histogramCalculator(formatedR1Residue1, formatedR1Residue2)
    residue1NOCoordinates = reshape(pagetranspose(mean(double(formatedR1Residue1(9:10, 7:9, :)))), 3, [])';
    residue2NOCoordinates = reshape(pagetranspose(mean(double(formatedR1Residue2(9:10, 7:9, :)))), 3, [])';
    pairDistanceBetweenNO1AndNO2 = pdist2(residue1NOCoordinates, residue2NOCoordinates);
    pairDistanceBetweenNO1AndNO2 = reshape(pairDistanceBetweenNO1AndNO2, [], 1);
    [simulatedDistanceDistributionY, simulatedDistanceDistributionX] = hist(pairDistanceBetweenNO1AndNO2, [10:90]);
    normalizedSimulatedDistanceDistributionY = simulatedDistanceDistributionY/max(simulatedDistanceDistributionY);
end