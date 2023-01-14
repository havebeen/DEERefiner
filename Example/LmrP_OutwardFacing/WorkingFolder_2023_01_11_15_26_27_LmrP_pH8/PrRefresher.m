function PrRefresher(app)
    PrNumber = length(app.DEERRefineRestraintsTable.Data(:, 1));
    app.DEERRefineRestraintsTable.Data(:, 1) = 1:PrNumber;
    app.numberOfTimesPrLoadingButtonHit = PrNumber;
end