function rowOfTableRemovingColorer(app)
    if app.rowOfTableToBeColored > size(app.DEERRefineRestraintsTable.Data, 1)
        return
    elseif app.rowOfTableToBeColored==0
        return
    else
        tableRemovingColorStyle = uistyle;
        tableRemovingColorStyle.BackgroundColor = "#EDB120";
        addStyle(app.DEERRefineRestraintsTable,tableRemovingColorStyle, 'row', app.rowOfTableToBeColored);
    end
end