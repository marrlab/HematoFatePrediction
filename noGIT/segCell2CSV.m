basepath = 'E:';
moviename = '140206PH8';
treeNr = 55;

outfile = 'F:\Data\Philipp\FatePrediction_2013\NatMethDemo\tracks_localXY.csv';
cellfiles = dir(sprintf('%s/%s/treemeasurements/tree_%.4i/track*.mat',basepath,moviename,treeNr));


fh = fopen(outfile,'w');
fprintf(fh,'pos;frame;x;y;cellnr\n');
for cf = cellfiles'
    c = load(sprintf('%s\\%s\\treemeasurements\\tree_%.4i\\%s',basepath,moviename,treeNr,cf.name),'segCell','celltrack');
    fn = c.segCell.featNames;
    assert(size(c.celltrack.positionIndex,2) == size(c.segCell.dataPoints,1))
    
    for i = 1:size(c.segCell.dataPoints,1)
        fprintf(fh,'%i;%i;%i;%i;%i;-1.0\n',c.celltrack.positionIndex(i),...
                                        c.segCell.dataPoints(i,strcmp(fn,'tp')),...
                                        c.celltrack(i,
                                        %c.segCell.dataPoints(i,strcmp(fn,'X')),...
                                        %c.segCell.dataPoints(i,strcmp(fn,'Y')),...
                                        c.segCell.nr);
        clear segCell celltrack
    end
    
end
fclose(fh);
