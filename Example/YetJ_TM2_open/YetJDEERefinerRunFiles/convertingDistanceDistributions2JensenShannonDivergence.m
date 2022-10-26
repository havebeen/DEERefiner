function JensenShannonDivergence = convertingDistanceDistributions2JensenShannonDivergence(distanceDistributionOne, distanceDistributionTwo)
    yShiftedDistanceDistributionOne = distanceDistributionOne+0.00000000001;
    yShiftedDistanceDistributionTwo = distanceDistributionTwo+0.00000000001;
    normalizedyShiftedDistanceDistributionOne = yShiftedDistanceDistributionOne/sum(yShiftedDistanceDistributionOne);
    normalizedyShiftedDistanceDistributionTwo = yShiftedDistanceDistributionOne/sum(yShiftedDistanceDistributionTwo);
    distanceDistributionsMean = (normalizedyShiftedDistanceDistributionOne+normalizedyShiftedDistanceDistributionTwo)/2;
    JensenShannonDivergenceFirstTerm = 1/2*sum(normalizedyShiftedDistanceDistributionOne.*log2(normalizedyShiftedDistanceDistributionOne./distanceDistributionsMean));
    JensenShannonDivergenceSecondTerm = 1/2*sum(normalizedyShiftedDistanceDistributionTwo.*log2(normalizedyShiftedDistanceDistributionTwo./distanceDistributionsMean));
    JensenShannonDivergence = JensenShannonDivergenceFirstTerm+JensenShannonDivergenceSecondTerm;
    
    
end