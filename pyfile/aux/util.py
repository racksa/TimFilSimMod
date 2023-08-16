import math

def spherical_to_cartesian(r, theta, phi):
    """
    Convert spherical coordinates to Cartesian coordinates.
    
    Args:
        r (float): Radius or distance from origin to point.
        theta (float): Polar angle in radians (0 <= theta <= pi).
        phi (float): Azimuthal angle in radians (0 <= phi <= 2*pi).
    
    Returns:
        tuple: Cartesian coordinates (x, y, z).
    """
    x = r * math.sin(theta) * math.cos(phi)
    y = r * math.sin(theta) * math.sin(phi)
    z = r * math.cos(theta)
    return x, y, z