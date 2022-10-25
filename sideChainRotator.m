function sideChainRotatedFormatedStructure = sideChainRotator(formatedStructure, rotatingAngle, rotatingResidueIndex)
    residueIndex = unique(double(formatedStructure(:, 7)));
    
    sideChainMatrix = strings(14, 16, length(residueIndex));
    sideChainFirstRotation = strings(14, 16, length(residueIndex));
    sideChainSecondRotation = strings(14, 16, length(residueIndex));
    sideChainThirdRotation = strings(14, 16, length(residueIndex));
    sideChainFourthRotation = strings(14, 16, length(residueIndex));
    sideChainFifthRotation = strings(14, 16, length(residueIndex));
    
    for i = 1:length(residueIndex)
        sideChainMatrix(1:length(formatedStructure(double(formatedStructure(:, 7)) == (residueIndex(i)), 1)), :, i) = ...
            (formatedStructure(double(formatedStructure(:, 7)) == (residueIndex(i)), :));
        if sum(residueIndex(i)==rotatingResidueIndex) > 0
            sideChainFirstRotation(:, :, i) = sideChainFirstBondRotator(sideChainMatrix(:, :, i), rotatingAngle);
            sideChainSecondRotation(:, :, i) = sideChainSecondBondRotator(sideChainFirstRotation(:, :, i), rotatingAngle);
            sideChainThirdRotation(:, :, i) = sideChainThirdBondRotator(sideChainSecondRotation(:, :, i), rotatingAngle);
            sideChainFourthRotation(:, :, i) = sideChainFourthBondRotator(sideChainThirdRotation(:, :, i), rotatingAngle);
            sideChainFifthRotation(:, :, i) = sideChainFifthBondRotator(sideChainFourthRotation(:, :, i), rotatingAngle);
        else
            sideChainFifthRotation(:, :, i) = sideChainMatrix(:, :, i);
        end
    end
    
    reshapedSideChainFifthRotation = reshape(pagetranspose(sideChainFifthRotation), 16, [])';
    sideChainRotatedFormatedStructure = reshapedSideChainFifthRotation(strlength(reshapedSideChainFifthRotation(:, 1)) ~= 0, :);
    sideChainRotatedFormatedStructure(:, 2) = string(1:length(sideChainRotatedFormatedStructure(:, 2)));
end