function clashesAutoInputer(app)
    residueIndex = unique(double(app.initialMutatedFormatedPDB(:, 7)));
    residueNumber = length(residueIndex);
    clashesNumber = round(residueNumber*2/100);
    app.MaximalclashesEditField.Value = clashesNumber;

end