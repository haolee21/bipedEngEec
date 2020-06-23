% input parameter style
% DH is a matrix
% DH = [alpha,a,d,theta,0,0]
% rad = 0 if the input angle will be in degree

function Robot = DH2Robot(DH,rad)

[m,n] = size(DH);

for i=1:m
    a = DH(i,2);
    d = DH(i,3);
    rot = DH(i,5);
    off = 0;
    al = 0;
    th = 0;
    
    if(rad == 0)
        al = DH(i,1)*pi/180;
        th = DH(i,4)*pi/180;
        if(rot == 0)
           off = DH(i,6)*pi/180;
        else
           off = DH(i,6);
        end
    else
        al = DH(i,1);
        th = DH(i,4);
        off = DH(i,6);
    end
    
    if(off == 0)
    L(i) = Link([th,d,a,al,rot],'modified');
    else
    L(i) = Link([th,d,a,al,rot,off],'modified');
    end
end

Robot = SerialLink(L);
end