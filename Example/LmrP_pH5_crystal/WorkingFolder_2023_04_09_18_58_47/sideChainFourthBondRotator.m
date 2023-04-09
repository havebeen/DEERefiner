function fourthBondRotatedSideChain = sideChainFourthBondRotator(formatedSideChain, rotatingAngle)
    residueType = formatedSideChain(1, 5);
    if residueType == "GLY" || residueType == "PRO" || residueType == "ALA" || ...
        residueType == "VAL" || residueType == "CYS" || residueType == "SER" || ...
        residueType == "THR" || residueType == "ILE" || residueType == "LEU" || ...
        residueType == "PHE" || residueType == "TYR" || residueType == "TRP" || ...
        residueType == "ASN" || residueType == "HIS" || residueType == "ASP" ||...
        residueType == "GLU" || residueType == "GLN" || residueType == "MET"
        fourthBondRotatedSideChain = formatedSideChain;
    elseif residueType == "ARG" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCDCoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CD ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCDCoordinates;
        ARGRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " NE ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CZ " |...
                     formatedSideChain(:, 3) == " NH1" |...
                     formatedSideChain(:, 3) == " NH2", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CZ " |...
                     formatedSideChain(:, 3) == " NH1" |...
                     formatedSideChain(:, 3) == " NH2", :)*ARGRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCDCoordinates;
        fourthBondRotatedSideChainTmp = formatedSideChain;
        fourthBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        fourthBondRotatedSideChain = fourthBondRotatedSideChainTmp; 
    elseif residueType == "LYS" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCDCoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CD ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCDCoordinates;
        LYSRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CE ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);  
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " NZ ", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " NZ ", :)*LYSRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCDCoordinates;
        fourthBondRotatedSideChainTmp = formatedSideChain;
        fourthBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        fourthBondRotatedSideChain = fourthBondRotatedSideChainTmp; 
    end
end