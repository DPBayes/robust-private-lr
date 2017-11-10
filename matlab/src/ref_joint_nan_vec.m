function [ref]=ref_joint_nan_vec(x,y)
    z=x+y;
    i=find(isnan(z)==0);
    x1=x(i);
    y1=y(i);
    ref.x=x1;
    ref.y=y1;
    
