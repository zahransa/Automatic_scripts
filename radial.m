function vector = radial(electrode1, electrode2, electrode3, dipodeCenter)
    
    x1 = electrode1(1);
    y1 = electrode1(2);
    z1 = electrode1(3);
    x2 = electrode2(1);
    y2 = electrode2(2);
    z2 = electrode2(3);
    x3 = electrode3(1);
    y3 = electrode3(2);
    z3 = electrode3(3);
    x4 = dipodeCenter(1);
    y4 = dipodeCenter(2);
    z4 = dipodeCenter(3);
    
    A = y1*(z2-z3)+y2*(z3-z1)+y3*(z1-z2);
    B = z1*(x2-x3)+z2*(x3-x1)+z3*(x1-x2);
    C = x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2);
    D = -x1*(y2*z3-y3*z2)-x2*(y3*z1-y1*z3)-x3*(y1*z2-y2*z1);
    
    A1 = A/D;
    B1 = B/D;
    C1 = C/D;
%     D1 = 1;
    
    mat1 = [x2-x1, y2-y1, z2-z1; ...
        x3-x1, y3-y1, z3-z1; ...
        A1, B1, C1];
    mat2 = [x4*(x2-x1)+y4*(y2-y1)+z4*(z2-z1);...
        x4*(x3-x1)+y4*(y3-y1)+z4*(z3-z1);...
        -1];
    vector = mat1\mat2;
    
    length = norm(vector);
    vector = vector/length;
        
end