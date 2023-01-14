function batchJob = runFilesClusterExecutor()
    runFiles = dir('*runFile_*');
    for execution = 1:length(runFiles)
        runFileNameTmp = erase(runFiles(execution).name, '.m');
        batchJob(execution) = batch(runFileNameTmp);
        pause(3);
    end
end