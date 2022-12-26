function [distanceDistributionFullPath, distanceDistributionFileName] = loadingDistanceDistribution
    [distanceDistributionFileNameTmp, distanceDistributionPathName] = uigetfile('*.dat;*.pr2');
    if distanceDistributionFileNameTmp==0
        distanceDistributionFullPath = "";
        distanceDistributionFileName = "";
      return
    end
    distanceDistributionFileName = string(distanceDistributionFileNameTmp);
    cd(distanceDistributionPathName)
    distanceDistributionFullPath = strcat(distanceDistributionPathName, distanceDistributionFileName);
end