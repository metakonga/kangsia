function plot_isph(p)
len = length(p);
% colors = size(len);
% ncmap=length(cmap);
% diff = pmax / ncmap;
% for i = 1:1:len
%     for j = 1:1:ncmap
%         vmin=diff*(j - 1);
%         vmax=diff*(j);
%         if(p(i,6) > vmin && p(i,6) <= vmax)
%             colors(i)=(j);
%             break;
%         end
%     end
% end
for i = 1:1:len
%     if(p(i,2) == 1)
%         plot(p(i,3), p(i,4), 'b', 'Marker', 'x', 'MarkerSize', 3);
%          hold on;
%         continue;
%     end
% %     c=cmap(colors(i),:);
%     if(p(i,3) > 10)
%         continue;
%     end
    if(p(i,7) > 0)
        plot(p(i,4), p(i,5), 'MarkerFaceColor', 'k', 'Marker', 'x', 'MarkerSize', 4);
        continue;
    end

    if(p(i,3) == 1)
        plot(p(i,4), p(i,5),'MarkerFaceColor', 'b', 'Marker', 'o', 'MarkerSize', 4);
    elseif(p(i,3) == 2)
        plot(p(i,4), p(i,5),'MarkerFaceColor', 'r', 'Marker', 'o', 'MarkerSize', 4);
    elseif(p(i,3) == 3)
        plot(p(i,4), p(i,5),'MarkerFaceColor', 'g', 'Marker', 'o', 'MarkerSize', 4);
    elseif(p(i,3) == 5)
        plot(p(i,4), p(i,5),'MarkerFaceColor', 'b', 'Marker', 'x', 'MarkerSize', 4);
    end
   
    
        
    hold on;
end