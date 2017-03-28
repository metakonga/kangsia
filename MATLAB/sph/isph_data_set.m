function [minval, maxval, tj, ti] = isph_data_set(data, t, cmap, dim, isvm)
nrow = size(data, 1);
id = zeros(nrow,1);
ty = zeros(nrow,1);
x = zeros(nrow,1);
y = zeros(nrow,1);
z = zeros(nrow,1);
v = zeros(nrow,1);
p = zeros(nrow,1);
% //out = zeros(nrow, 8);
for i = 1:nrow
    id(i) = data(i, 2);
    ty(i) = data(i, 3);
    x(i) = data(i, 4);
    y(i) = data(i, 5);
    z(i) = data(i, 6);
    if(isvm == 0)
        v(i) = data(i, 7); 
    else
        v(i) = sqrt(data(i, 7) * data(i, 7) + data(i, 8) * data(i, 8) + data(i, 9) * data(i, 9));
    end
    p(i) = data(i, 10);
%     out(i, 8) = data(i, 11);
end
maxval = 0; 
ti = 0;
if(t == 'v')
    [maxval, ti] = max(v);
    [minval, tj] = min(v);
    plot_isph_particles(minval, maxval, ty, x, y, z, v, cmap, dim);
elseif(t == 'p')
    [maxval, ti] = max(p);
    [minval, tj] = min(p);
    plot_isph_particles(minval, maxval, ty, x, y, z, p, cmap, dim);
end

% for i = 1:nrow
%     if(out(i,2) == 1)
%         plot(out(i,3), out(i,4),'MarkerFaceColor', 'b', 'Marker', 'o', 'MarkerSize', 4);
%     elseif(out(i,2) == 2)
%         plot(out(i,3), out(i,4),'MarkerFaceColor', 'r', 'Marker', 'o', 'MarkerSize', 4);
%     elseif(out(i,2) == 3)
%         plot(out(i,3), out(i,4),'MarkerFaceColor', 'g', 'Marker', 'o', 'MarkerSize', 4);
%     elseif(out(i,2) == 5)
%         plot(out(i,3), out(i,4),'MarkerFaceColor', 'b', 'Marker', 'x', 'MarkerSize', 4);
%     end
%     hold on;
% end