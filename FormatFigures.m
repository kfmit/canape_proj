function FormatFigures(NameFig, save_mat_fig)


    
title('')

box on

h = get(0,'children');
scrsz = get(0,'ScreenSize');
set(h,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)])  

set(gcf,'color','w'); 

print(gcf,NameFig,'-dpng')
if save_mat_fig
    saveas(gcf,NameFig,'fig')
end
    
close(gcf)

end
