function sortingDEERRefineRestraintsTable(app)
    app.DEERRefineRestraintsTable.Data(:, 1) = string(1:size(app.DEERRefineRestraintsTable.Data, 1))';
end