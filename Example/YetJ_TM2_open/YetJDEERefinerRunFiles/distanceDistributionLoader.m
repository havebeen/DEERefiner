function DEERAnalysisData = distanceDistributionLoader(app)
    PrFileType = PrFileTypeReader(app.distanceDistributionFullPath(app.numberOfTimesPrLoadingButtonHit));
    if PrFileType== "DAPr"
        DEERAnalysisData = load(app.distanceDistributionFullPath(app.numberOfTimesPrLoadingButtonHit));
    elseif PrFileType== "pr2"
        fid = fopen(app.distanceDistributionFullPath(app.numberOfTimesPrLoadingButtonHit), 'r');
        sf_cell = textscan(fid, '%f', 'CommentStyle', '#');
        fclose all;
        sf_double = sf_cell{1};
        sf_double(1:2) = [];
        reshaped_sf_double = reshape(sf_double, [], 1024)';
        DEERAnalysisData = reshaped_sf_double(:, 1:2);
    end
    
end