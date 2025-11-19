function v_global = rotate_vector_fault2global(v_local, strike_deg, dip_deg)

 % v_local: N×3 matrix (each row is [strike_comp, dip_comp, normal_comp])
 % strike_deg, dip_deg: scalars (same orientation for all vectors)

    % Convert angles
    strike = deg2rad(strike_deg);
    dip = deg2rad(dip_deg);

    % Compute local basis vectors in global Cartesian (Z up, Y north, X east)
    e_strike = [sin(strike), cos(strike), 0];  % horizontal, strike direction
    e_dip    = [cos(strike).*cos(dip),  -sin(strike).*cos(dip), -sin(dip)];  
    e_normal = cross( e_dip,e_strike);  % right-handed system

    % Build rotation matrix: columns are local basis vectors
    R = [ e_dip(:)/norm(e_dip), e_strike(:)/norm(e_strike),e_normal(:)/norm(e_normal)];

    % Multiply all vectors (matrix product)
    v_global = v_local * R';  % N×3 * 3×3' = N×3
end