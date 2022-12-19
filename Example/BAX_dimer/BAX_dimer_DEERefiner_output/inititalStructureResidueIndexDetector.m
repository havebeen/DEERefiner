function residueIndex = inititalStructureResidueIndexDetector(app)
    residueIndex = unique(double(app.initialMutatedFormatedPDB(:, 7)));
    app.FlexibleregionsEditField.FontColor = [0 0 0];
end