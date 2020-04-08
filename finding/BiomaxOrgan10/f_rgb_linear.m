function rgb = f_rgb_linear(XYZ)
    % constants
    m=[3.2410 -1.5374 -0.4986; -0.9692 1.8760 0.0416; 0.0556 -0.2040 1.0570];
    
    % linearize
    rgb = m*XYZ';
end

