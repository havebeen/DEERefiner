function loopRegionString = loopRegionDeterminer(app)
    loopRegion = unique(double(app.initialMutatedFormatedPDB(:, 7)));
    loopRegionStringStart = loopRegion(1);
    loopRegionStringEnd = loopRegion(end);
    loopRegionString = strcat(num2str(loopRegionStringStart), ":", num2str(loopRegionStringEnd-1));
end