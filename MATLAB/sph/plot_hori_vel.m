function plot_hori_vel(p, min_r, max_r)
nrow = size(p, 1);
px = zeros(nrow, 1);
v = zeros(nrow,1);
cnt = 1;
for i = 1:nrow
    if(p(i,4) > min_r && p(i,4) < max_r)
        px(cnt) = p(i, 5);
        v(cnt) = p(i,7);
        cnt = cnt + 1;
         plot(p(i,7), p(i,5), 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'Marker', 'o', 'MarkerSize', 6);
         hold on;
    end
end

