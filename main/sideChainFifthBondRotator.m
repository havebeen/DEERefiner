function fifthBondRotatedSideChain = sideChainFifthBondRotator(formatedSideChain, rotatingAngle)
    residueType = formatedSideChain(1, 5);
    if residueType == "GLY" || residueType == "PRO" || residueType == "ALA" || ...
        residueType == "VAL" || residueType == "CYS" || residueType == "SER" || ...
        residueType == "THR" || residueType == "ILE" || residueType == "LEU" || ...
        residueType == "PHE" || residueType == "TYR" || residueType == "TRP" || ...
        residueType == "ASN" || residueType == "HIS" || residueType == "ASP" ||...
        residueType == "GLU" || residueType == "GLN" || residueType == "MET" ||...
        residueType == "LYS"
        fifthBondRotatedSideChain = formatedSideChain;
    elseif residueType == "ARG" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetNECoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " NE ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetNECoordinates;
        ARGRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CZ ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " NH1" |...
                     formatedSideChain(:, 3) == " NH2", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " NH1" |...
                     formatedSideChain(:, 3) == " NH2", :)*ARGRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetNECoordinates;
        fifthBondRotatedSideChainTmp = formatedSideChain;
        fifthBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        fifthBondRotatedSideChain = fifthBondRotatedSideChainTmp; 
    end
end