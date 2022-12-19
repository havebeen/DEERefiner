function distanceDistribution = distanceDistributionFinalLoader(distanceDistributionFileName)
    PrFileType = PrFileTypeReader(distanceDistributionFileName);
    if PrFileType== "DAPr"
        distanceDistributionRaw = load(distanceDistributionFileName);
        distanceDistributionTmp(:, 1) = distanceDistributionRaw(:, 1)*10;
        distanceDistributionTmp(:, 2) = distanceDistributionRaw(:, 2)/sum(distanceDistributionRaw(:, 2));
        distanceDistributionX = 10:90;
        unnormalizedDistanceDistribution = interp1(distanceDistributionTmp(:, 1), distanceDistributionTmp(:, 2), distanceDistributionX);
        unnormalizedDistanceDistribution(isnan(unnormalizedDistanceDistribution)) = 0;
        distanceDistributionY = unnormalizedDistanceDistribution/sum(unnormalizedDistanceDistribution);
        distanceDistribution = [distanceDistributionX' distanceDistributionY'];
    elseif PrFileType== "pr2"
        fid = fopen(distanceDistributionFileName, 'r');
        sf_cell = textscan(fid, '%f', 'CommentStyle', '#');
        fclose all;
        sf_double = sf_cell{1};
        sf_double(1:2) = [];
        reshaped_sf_double = reshape(sf_double, [], 1024)';
        distanceDistributionTmp(:, 1) = reshaped_sf_double(:, 1)*10;
        distanceDistributionTmp(:, 2) = reshaped_sf_double(:, 2)/sum(reshaped_sf_double(:, 2));
        distanceDistributionX = 10:90;
        unnormalizedDistanceDistribution = interp1(distanceDistributionTmp(:, 1), distanceDistributionTmp(:, 2), distanceDistributionX);
        unnormalizedDistanceDistribution(isnan(unnormalizedDistanceDistribution)) = 0;
        distanceDistributionY = unnormalizedDistanceDistribution/sum(unnormalizedDistanceDistribution);
        distanceDistribution = [distanceDistributionX' distanceDistributionY'];
    end
end