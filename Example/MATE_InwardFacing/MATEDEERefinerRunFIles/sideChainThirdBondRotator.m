function thirdBondRotatedSideChain = sideChainThirdBondRotator(formatedSideChain, rotatingAngle)
    residueType = formatedSideChain(1, 5);
    if residueType == "GLY" || residueType == "PRO" || residueType == "ALA" || ...
        residueType == "VAL" || residueType == "CYS" || residueType == "SER" || ...
        residueType == "THR" || residueType == "ILE" || residueType == "LEU" || ...
        residueType == "PHE" || residueType == "TYR" || residueType == "TRP" || ...
        residueType == "ASN" || residueType == "HIS" || residueType == "ASP"
        thirdBondRotatedSideChain = formatedSideChain;
    elseif residueType == "ARG" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCGCoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CG ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCGCoordinates;
        ARGRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CD ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " NE " |...
                     formatedSideChain(:, 3) == " CZ " |...
                     formatedSideChain(:, 3) == " NH1" |...
                     formatedSideChain(:, 3) == " NH2", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " NE " |...
                     formatedSideChain(:, 3) == " CZ " |...
                     formatedSideChain(:, 3) == " NH1" |...
                     formatedSideChain(:, 3) == " NH2", :)*ARGRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCGCoordinates;
        thirdBondRotatedSideChainTmp = formatedSideChain;
        thirdBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        thirdBondRotatedSideChain = thirdBondRotatedSideChainTmp; 
    elseif residueType == "LYS" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCGCoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CG ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCGCoordinates;
        LYSRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CD ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CE " |...
                     formatedSideChain(:, 3) == " NZ ", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CE " |...
                     formatedSideChain(:, 3) == " NZ ", :)*LYSRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCGCoordinates;
        thirdBondRotatedSideChainTmp = formatedSideChain;
        thirdBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        thirdBondRotatedSideChain = thirdBondRotatedSideChainTmp;
    elseif residueType == "GLU" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCGCoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CG ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCGCoordinates;
        GLURotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CD ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " OE1" |...
                     formatedSideChain(:, 3) == " OE2", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " OE1" |...
                     formatedSideChain(:, 3) == " OE2", :)*GLURotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCGCoordinates;
        thirdBondRotatedSideChainTmp = formatedSideChain;
        thirdBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        thirdBondRotatedSideChain = thirdBondRotatedSideChainTmp;
    elseif residueType == "GLN" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCGCoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CG ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCGCoordinates;
        GLNRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CD ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " OE1" |...
                     formatedSideChain(:, 3) == " NE2", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " OE1" |...
                     formatedSideChain(:, 3) == " NE2", :)*GLNRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCGCoordinates;
        thirdBondRotatedSideChainTmp = formatedSideChain;
        thirdBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        thirdBondRotatedSideChain = thirdBondRotatedSideChainTmp;
    elseif residueType == "MET" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCGCoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CG ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCGCoordinates;
        METRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " SD ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CE ", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CE ", :)*METRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCGCoordinates;
        thirdBondRotatedSideChainTmp = formatedSideChain;
        thirdBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        thirdBondRotatedSideChain = thirdBondRotatedSideChainTmp;

    end
end