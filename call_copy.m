function call_copy()
% clicking on axes within a figure, opens a new figure with only the
% clicked axes in it.
% Usage:
% ff = findobj('type','axes');
% set(ff,'ButtonDownFcn','call_copy');

% (c) Tonio

copyobj(gca,figure);
set(gca,'Position',[0.13 0.11 0.775 0.815]);
set(gca,'ButtonDownFcn','');
h_child=get(gca,'Children');
set(h_child,'ButtonDownFcn','');
if strcmp(get(h_child(1),'Type'),'surface') || strcmp(get(h_child(end),'Type'),'image')
    colorbar;
end;
set(gca, 'XTickMode','auto');
set(gca, 'YTickMode','auto');
