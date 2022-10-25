function minimalBackboneNonBondedDistance = formatedBackbone2MinimalBackboneNonBondedDistance(formatedBackbone)
    backboneCoordinates = double(formatedBackbone(:, 9:11));
    backboneCoordinatesPairDistance = pdist2(backboneCoordinates, backboneCoordinates);
    bondedFlag = eye(length(backboneCoordinates(:, 1)));
    residueIndex = unique(double(formatedBackbone(:, 7)));
    for loopResidueIndex = 1:length(residueIndex)
        bondedFlag(formatedBackbone(:, 3) == " N  " & double(formatedBackbone(:, 7)) == residueIndex(loopResidueIndex), formatedBackbone(:, 3) == " CA " & double(formatedBackbone(:, 7)) == residueIndex(loopResidueIndex)) = 1;
        bondedFlag(formatedBackbone(:, 3) == " CA " & double(formatedBackbone(:, 7)) == residueIndex(loopResidueIndex), formatedBackbone(:, 3) == " C  " & double(formatedBackbone(:, 7)) == residueIndex(loopResidueIndex)) = 1;
        bondedFlag(formatedBackbone(:, 3) == " C  " & double(formatedBackbone(:, 7)) == residueIndex(loopResidueIndex), formatedBackbone(:, 3) == " O  " & double(formatedBackbone(:, 7)) == residueIndex(loopResidueIndex)) = 1;
        bondedFlag(formatedBackbone(:, 3) == " N  " & double(formatedBackbone(:, 7)) == residueIndex(loopResidueIndex), formatedBackbone(:, 3) == " C  " & double(formatedBackbone(:, 7)) == residueIndex(loopResidueIndex)) = 1;
        bondedFlag(formatedBackbone(:, 3) == " CA " & double(formatedBackbone(:, 7)) == residueIndex(loopResidueIndex), formatedBackbone(:, 3) == " O  " & double(formatedBackbone(:, 7)) == residueIndex(loopResidueIndex)) = 1;
    end
    for loopResidueIndex = 1:length(residueIndex)-1
        bondedFlag(formatedBackbone(:, 3) == " C  " & double(formatedBackbone(:, 7)) == residueIndex(loopResidueIndex), formatedBackbone(:, 3) == " N  " & double(formatedBackbone(:, 7)) == residueIndex(loopResidueIndex+1)) = 1;
        bondedFlag(formatedBackbone(:, 3) == " O  " & double(formatedBackbone(:, 7)) == residueIndex(loopResidueIndex), formatedBackbone(:, 3) == " N  " & double(formatedBackbone(:, 7)) == residueIndex(loopResidueIndex+1)) = 1;
        bondedFlag(formatedBackbone(:, 3) == " CA " & double(formatedBackbone(:, 7)) == residueIndex(loopResidueIndex), formatedBackbone(:, 3) == " N  " & double(formatedBackbone(:, 7)) == residueIndex(loopResidueIndex+1)) = 1;
        bondedFlag(formatedBackbone(:, 3) == " C  " & double(formatedBackbone(:, 7)) == residueIndex(loopResidueIndex), formatedBackbone(:, 3) == " CA " & double(formatedBackbone(:, 7)) == residueIndex(loopResidueIndex+1)) = 1;
    end
    bondedFlag = bondedFlag+triu(bondedFlag, 1)';
    finalBackboneCoordinatesPairDistance = backboneCoordinatesPairDistance+5000*bondedFlag;
    minimalBackboneNonBondedDistance = min(min(finalBackboneCoordinatesPairDistance));
end