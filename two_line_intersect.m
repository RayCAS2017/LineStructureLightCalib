function intersect = two_line_intersect(line1, line2)
% line1 = [a,c] = ax + c
% line2 = [b,d] = bx + d

intersect_x = (line2(2) - line1(2))/(line1(1)-line2(1));
intersect_y = line1(1)*intersect_x + line1(2);
intersect=[intersect_x,intersect_y];