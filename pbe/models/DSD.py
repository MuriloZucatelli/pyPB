    """Module for droplet size distribuition functions
    """

#importing dependences
from numpy import exp, pi

def A0(self, v):
    """_summary_

    Args:
        v (_type_): _description_

    Returns:
        _type_: _description_
    """
    return \
        self.phi / (self.v0 * self.sigma0 * sqrt(2 * pi)) * \
        exp(-(v - self.v0)**2 / (2 * self.sigma0**2))