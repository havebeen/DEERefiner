function PrFileType = PrFileTypeReader(PrFileName)
    try
        Pr = load(PrFileName);
        PrFileType = "DAPr";
    catch
        PrFileType = "pr2";
    end
end