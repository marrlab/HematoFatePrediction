moviepath = 'E:\';
moviename = '140206PH8';
movie_identifier = 'experiment3';
treenr = 55;
branchnrs = [164,193];

load([moviepath 'resultfiles\results' moviename])
outpath = 'F:\Data\Philipp\FatePrediction_2013\NatMethDemo\';
mkdir(outpath);
bgoutpath = [outpath movie_identifier '\backgrounds\'];
mkdir(bgoutpath);

%%
allcells = [];
for br = 1:numel(branchnrs)
    allcells = [allcells rootpath(branchnrs(br))];
end
allcells = unique(allcells);

posis = unique(results{treenr}.nonFluor.positionIndex(ismember(results{treenr}.nonFluor.cellNr,allcells)));
mintp = min(results{treenr}.nonFluor.timepoint(ismember(results{treenr}.nonFluor.cellNr,allcells)));
maxtp = max(results{treenr}.nonFluor.timepoint(ismember(results{treenr}.nonFluor.cellNr,allcells)));


%%
track = results{treenr}.nonFluor;
fieldz  = fields(track);
keepfieldz = {'cellNr','timepoint','absoluteX','absoluteY','positionIndex','X',...
    'Y','wavelength_3','wavelength_2','wavelength_1','absoluteTime','filename','invtype_alt',...
    'invgen_alt','type'};

track = rmfield(track,fieldz(~ismember(fieldz,keepfieldz)));

fieldz = fields(track);
for f = fieldz'
    if strcmp(f{:},'cellNr')
        continue
    end
        track.(f{:})(~ismember(track.cellNr,allcells)) = [];
end
track.cellNr(~ismember(track.cellNr,allcells)) = [];

for t = 1:numel(track.filename)
   [~,tokens] = regexp(track.filename(t),'^.*/(.*)_p(\d+)_t(\d+)_z001_w(\d+).png$','match','tokens');
   newname = sprintf('pos%s_t%s_wl%s.png',tokens{1}{1}{2},tokens{1}{1}{3},tokens{1}{1}{4});
   track.filename{t} = newname;
end
save([outpath 'track.mat'],'track')
%%
mkdir([outpath '\' movie_identifier])

for p=posis
    mkdir([outpath '\' movie_identifier '\position' num2str(p)])
    mkdir([bgoutpath '\position' num2str(p)])
    projbgcopied = 0;
    for t = mintp:maxtp
        for wl = [0,1,2,3];
%             if exist(sprintf('%s%s/%s_p%.4i/%s_p%.4i_t%.5i_z001_w%.2i.png', moviepath, moviename, moviename, p, moviename, p, t, wl),'file')
%                 copyfile(sprintf('%s%s/%s_p%.4i/%s_p%.4i_t%.5i_z001_w%.2i.png', moviepath, moviename, moviename, p, moviename, p, t, wl),...
%                     sprintf('%s/%s/position%i/%s_pos%.4i_t%.5i_wl%.2i.png',outpath,movie_identifier,p,movie_identifier,p,t,wl));
%             end
           if exist(sprintf('%s%s/background/%s_p%.4i/%s_p%.4i_t%.5i_z001_w%.2i.png', moviepath, moviename, moviename, p, moviename, p, t, wl),'file')
                 copyfile(sprintf('%s%s/background/%s_p%.4i/%s_p%.4i_t%.5i_z001_w%.2i.png', moviepath, moviename, moviename, p, moviename, p, t, wl),...
                     sprintf('%s/position%i/%s_pos%.4i_t%.5i_wl%.2i.png',bgoutpath,p,movie_identifier,p,t,wl));
           end
            if exist(sprintf('%s%s/background_projected/%s_p%.4i/%s_p%.4i_w00_projbg.png', moviepath, moviename, moviename, p, moviename, p, t, wl),'file') &&...
                    ~projbgcopied
                 copyfile(sprintf('%s%s/background_projected/%s_p%.4i/%s_p%.4i_w00_projbg.png', moviepath, moviename, moviename, p, moviename, p, t, wl),...
                     sprintf('%s/position%i/%s_pos%.4i_wl%.2i_projbg.png',bgoutpath,p,movie_identifier,p,t,wl));
                 projbgcopied = 1;
            end
        end
    end
end