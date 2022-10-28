function currentStatusUpdator(currentStatusFileName, currentStatusFlag)
    currentStatusFileIndex = fopen(currentStatusFileName, 'w');
    fprintf(currentStatusFileIndex, '%d', currentStatusFlag);
    fclose(currentStatusFileIndex);
end