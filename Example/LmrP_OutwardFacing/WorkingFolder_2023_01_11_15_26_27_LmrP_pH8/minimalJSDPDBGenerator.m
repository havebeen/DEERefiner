function minimalJSDIndex = minimalJSDPDBGenerator(FINALPDBENSEMBLE, JensenShannonDivergence, timeStamp)
    FINALPDBFileName = strcat("FINALPDB_",  timeStamp, ".pdb");
    [~, minimalJSDIndex] = min(JensenShannonDivergence);
    pdbSaver(FINALPDBFileName, FINALPDBENSEMBLE(:, :, minimalJSDIndex));
end