function plot_isph_particles(pmin, pmax, id, x, y, z, v, cmap, dim)
len = length(v);
colors = zeros(len,3);
ncmap=length(cmap);
diff = (pmax - pmin) / ncmap;
for i = 1:1:len
    if(v(i) >= -0.07 && v(i) < 0.02)
        colors(i, :) = [0 0.686274509 1];
%     elseif(v(i) >= -0.15 && v(i) < -0.07)
%         colors(i, :) = [1 0.4 0];
%     elseif(v(i) >= -0.07 && v(i) < 0.02)
%         colors(i, :) = [1 0.6 0];
%     elseif(v(i) >= 0.02 && v(i) < 0.10)
%         colors(i, :) = [1 1 0];
%     elseif(v(i) >= 0.1 && v(i) < 0.18)
%         colors(i, :) = [0.6 1 0];
%     elseif(v(i) >= 0.18 && v(i) < 0.26)
%         colors(i, :) = [0.341176480054855 1 0];
%     elseif(v(i) >= 0.26 && v(i) < 0.34)
%         colors(i, :) = [0 1 0];
%     elseif(v(i) >= 0.34 && v(i) < 0.43)
%         colors(i, :) = [0 1 0.4];
%     elseif(v(i) >= 0.43 && v(i) < 0.51)
%         colors(i, :) = [0 1 0.6];
%     elseif(v(i) >= 0.51 && v(i) < 0.59)
%         colors(i, :) = [0 1 0.8];
%     elseif(v(i) >= 0.59 && v(i) < 0.67)
%         colors(i, :) = [0 1 1];
%     elseif(v(i) >= 0.67 && v(i) < 0.76)
%         colors(i, :) = [0 0.8 1];
%     elseif(v(i) >= 0.76 && v(i) < 0.84)
%         colors(i, :) = [0 0.6 1];
%     elseif(v(i) >= 0.84 && v(i) < 0.92)
%         colors(i, :) = [0 0.4 1];
%     elseif(v(i) >= 0.92)
%         colors(i, :) = [0 0 1];
    end
%     for j = 1:1:ncmap
%         vmin=pmin + diff*(j - 1);
%         vmax=pmin + diff*(j);
%         if(v(i) >= vmin && v(i) <= vmax)
%             colors(i)=(j);
%             break;
%         end
%     end
end

if (dim == 2)
    for i = 1:1:len
    if(id(i) == 1)
%         c=cmap(colors(i),:);
        c=colors(i,:);
        plot(x(i), y(i), 'MarkerEdgeColor', c, 'MarkerFaceColor', c, 'Marker', 'o', 'MarkerSize', 6);
%     elseif(id(i) == 2)
%         plot(x(i), y(i),'MarkerFaceColor', 'r', 'Marker', 'o', 'MarkerSize', 4);
%     elseif(id(i) == 3)
%         plot(x(i), y(i),'MarkerFaceColor', 'g', 'Marker', 'o', 'MarkerSize', 4);
%     elseif(id(i) == 5)
%         plot(x(i), y(i),'MarkerFaceColor', 'b', 'Marker', 'x', 'MarkerSize', 4);
    end
    hold on;
    end
elseif(dim == 3)
    for i = 1:1:len
    if(id(i) == 1)
        c=cmap(colors(i),:);
        plot3(x(i), y(i), z(i), 'MarkerEdgeColor', c, 'MarkerFaceColor', c, 'Marker', 'o', 'MarkerSize', 6);
%     elseif(id(i) == 2)
%         plot(x(i), y(i), z(i),'MarkerFaceColor', 'r', 'Marker', 'o', 'MarkerSize', 4);
%     elseif(id(i) == 3)
%         plot(x(i), y(i), z(i),'MarkerFaceColor', 'g', 'Marker', 'o', 'MarkerSize', 4);
%     elseif(id(i) == 5)
%         plot(x(i), y(i), z(i),'MarkerFaceColor', 'b', 'Marker', 'x', 'MarkerSize', 4);
    end
    hold on;
    end
end


