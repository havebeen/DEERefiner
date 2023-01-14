function RMSEDatSaver(datFileName, currentStep, currentRMSE)
    fileIndex = fopen(datFileName, 'a');
    fprintf(fileIndex, "%d\t%f\n", currentStep, currentRMSE);
    fclose(fileIndex);
    

end