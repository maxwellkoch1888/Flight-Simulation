import numpy as np 

def quat_to_euler(quat: np.array) -> np.array:
    tol = 1e-12

    # EXTRACT e0, ex, ey, AND ez FROM quat
    e0 = quat(0)
    ex = quat(1)
    ey = quat(2)
    ez = quat(3)

    # CALCULATE RESULT FOR GIMBLE LOCK
    res = e0*ey - ex*ez

    # CALCULATE THE EQUIVALENT EULER ANGLES
    if abs(res - 0.5) < tol:
        azimuth_angle = 0.0
        bank_angle = 2*np.asin(ex/np.cos(np.pi/4.0) + azimuth_angle)
        pitch_angle = np.pi/2.0
    elif abs(res + 0.5) < tol: 
        azimuth_angle = 0.0
        bank_angle = 2.0*np.asin(ex/np.cos(np.pi/4.0) - azimuth_angle)
        pitch_angle = -np.pi/2
    else:
        bank_angle = np.atan2(2.0*(e0*ex + ey*ez), (e0**2 + ez**2 - ex**2 - ey**2))
        pitch_angle = np.asin(2.0*(e0*ey - ex*ez))
        azimuth_angle = np.atan2(2.0*(e0*ez + ex*ey), (e0**2 + ex**2 - ey**2 - ez**2))

    # RETURN THE EULER ANGLE
    euler_angles = np.array([bank_angle, pitch_angle, azimuth_angle])
    return euler_angles

def eul_to_quat(euler:np.array) -> np.array:
    # EXTRACT THE PITCH, BANK, AND AZIMUTH ANGLES
    pitch_angle   = euler[1]
    bank_angle    = euler[0]
    azimuth_angle = euler[2]

    # BUILD THE SIN AND COSINE VARIABLES
    cb = np.cos(bank_angle/2.0)
    sb = np.sin(bank_angle/2.0)
    cp = np.cos(pitch_angle/2.0)
    sp = np.sin(pitch_angle/2.0)
    ca = np.cos(azimuth_angle/2.0)
    sa = np.sin(azimuth_angle/2.0)

    # CALCULATE e0, ex, ey, ez
    e0 = cb*cp*ca + sb*sp*sa
    ex = sb*cp*ca - cb*sp*sa
    ey = cb*sp*ca + sb*cp*sa
    ez = cb*cp*sa - sb*sp*ca

    # BUILD THE QUATERNION
    quat = np.array([e0, ex, ey, ez])          
    return quat
     
def quat_mult(quat_a: np.array, quat_b: np.array)-> np.array:
    # EXTRACT THE VALUES FROM THE TWO QUATERNIONS
    a0 = quat_a[1]
    ax = quat_a[2]
    ay = quat_a[3]
    az = quat_a[4]

    b0 = quat_b[1]
    bx = quat_b[2]
    by = quat_b[3]
    bz = quat_b[4]

    # BUILD THE 0, X, Y, AND Z COMPONENTS OF THE NEW QUATERNION
    quat_c = np.zeros(4,1)
    quat_c[1] = a0*b0 - ax*bx - ay*by - az*bz
    quat_c[2] = a0*bx + ax*b0 + ay*bz - az*by
    quat_c[3] = a0*by - ax*bz + ay*b0 + az*bx
    quat_c[4] = a0*bz + ax*by - ay*bx + az*b0

    return quat_c

def quat_norm(quat: np.array) -> np.array:
    # FIND QUATERNION MAGNITUDE
    quat_mag = (quat(1)**2 + quat(2)**2 + quat(3)**2 + quat(4)**2)**0.5

    # NORMALIZE EACH COMPONENT OF THE QUATERNION
    quat[1] = quat[1] / quat_mag
    quat[2] = quat[2] / quat_mag
    quat[3] = quat[3] / quat_mag
    quat[4] = quat[4] / quat_mag
    
    return quat

def quat_dependent_to_base(dependent_vec: np.array, quat: np.array) -> np.array:
    # BUILD VARIABLES FROM EQUATION 1.5.4 
    vec_quat = np.array([0.0, dependent_vec(1), dependent_vec(2), dependent_vec(3)])
    quat_conj = np.array([quat(1), -quat(2), -quat(3), -quat(4)])
    
    # FIRST ROTATION FROM 1.5.4
    temp_quat = quat_mult(vec_quat, quat_conj)

    # SECOND ROTATION FROM 1.5.4
    quat_sol = quat_mult(quat, temp_quat)
    base_vec = np.array([quat_sol(2), quat_sol(3), quat_sol(4)])    
    return base_vec

def quat_base_to_dependent(base_vec:np.array, quat: np.array) -> np.array:
    vec_quat = np.array([0.0, base_vec(1), base_vec(2), base_vec(3)])
    quat_conj = np.array([quat(1), -quat(2), -quat(3), -quat(4)])

    # FIRST ROTATION FROM 1.5.4
    temp_quat = quat_mult(vec_quat, quat)

    # SECOND ROTATION FROM 1.5.4
    quat_sol = quat_mult(quat_conj, temp_quat)
    dependent_vec = np.array([quat_sol(2), quat_sol(3), quat_sol(4)])    
    return dependent_vec