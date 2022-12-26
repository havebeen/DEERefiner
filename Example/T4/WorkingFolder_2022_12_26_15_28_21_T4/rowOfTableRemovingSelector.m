function rowOfTableToBeColored = rowOfTableRemovingSelector(app)
    numberOfRowToBeRemoved = string(app.SelectingRemovingSpinPairSpinner.Value);
    [rowOfTableToBeColored, ~] = find(app.DEERRefineRestraintsTable.Data(:, 1)==numberOfRowToBeRemoved);
end