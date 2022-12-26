function mutatedPDBfile = DEERefineInitialSturctureMutator(app)
    initialFullFormatedPDB = pdbLoader(app.initialStructureFullPath);
    mutatedPDBfile = brokenFormatedPDBFixer(initialFullFormatedPDB);
end