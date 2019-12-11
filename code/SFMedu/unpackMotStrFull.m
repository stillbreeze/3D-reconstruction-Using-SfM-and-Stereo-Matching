function [Mot,Str,fx,fy,px,py] = unpackMotStrFull(nCam,vec)

cut = 3*2*nCam;
fx   = vec(1);
fy = vec(2);
px = vec(3);
py = vec(4);
Mot = reshape(vec(5:(cut + 4)),3,2,[]);
Str = reshape(vec(cut+5:end),3,[]);

end