function [ ] = gcfsavepdf( filename )
% this function saves the current figure to a pdf, with prescribed input filename
% change screen resolution from default 72 ppi to 109 ppi if using a mac
% (i.e., Thunderbolt display, MacBookPro)
if ismac
    set(0,'ScreenPixelsPerInch',109);
end

set(gcf,'Units','Inches');
figpos = get(gcf, 'Position');
paperWidth = figpos(3);
paperHeight = figpos(4);
set(gcf, 'paperunits', 'Inches');
set(gcf, 'papersize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0 0 paperWidth paperHeight]);

print('-dpdf',filename);

end

