function datSaver(datFileName, currentStep, currentRMSE, currentClashes)
    fileIndex = fopen(datFileName, 'a');
    fprintf(fileIndex, "%d\t%f\t%d\n", currentStep, currentRMSE, currentClashes);
    fclose(fileIndex);
    

end