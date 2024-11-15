function plot_ellipse(a,b,cx,cy,angle,color)
%a: width in pixels
%b: height in pixels
%cx: horizontal center
%cy: vertical center
%angle: orientation ellipse in radians (ccw from +x axis)
%color: color code (e.g., 'r' or [0.4 0.5 0.1])

angle = 2*pi - angle;
r=0:0.1:2*pi+0.1;
p=[(a*cos(r))' (b*sin(r))'];

alpha=[cos(angle) -sin(angle)
       sin(angle) cos(angle)];
   
 p1=p*alpha;
 
patch(cx+p1(:,1),cy+p1(:,2),color,'EdgeColor',color);

 
   
   