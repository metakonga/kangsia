function plot_pos(p)
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
    if(i > 5002)
        plot(p(i,1), p(i,3),'MarkerFaceColor', 'b', 'Marker', 'o', 'MarkerSize', 2);
    else
        plot(p(i,1), p(i,3),'MarkerFaceColor', 'r', 'Marker', 'o', 'MarkerSize', 4);
    end
%     //elseif(p(i,1) == 2)
%         plot(p(i,3), p(i,4),'MarkerFaceColor', 'r', 'Marker', 'o', 'MarkerSize', 4);
%    // elseif(p(i,1) == 3)
%         plot(p(i,3), p(i,4),'MarkerFaceColor', 'g', 'Marker', 'o', 'MarkerSize', 4);
%     end
   
        
    hold on;
end