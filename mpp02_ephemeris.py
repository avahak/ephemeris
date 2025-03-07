"""
Calculates position and velocity of the Moon w.r.t. the Earth using ELP/MPP02.
Can use LLR or DE405 fit data and their truncated versions.

ELP/MMP02 article: Lunar solution ELP, version ELP/MPP02, Jean Chapront and Gerard Francou. 
"""
import numpy as np

import tools.misc as misc

class MPP02Ephemeris:
    """
    Computes position and velocity of the Moon using ELP/MPP02.
    Units: km for position, km/Julian day for velocity.
    """
    ECL_TO_EQU = misc.rotation_matrix(0, 84381.448 * misc.ARCSEC)  # angle ~23.44 deg

    def __init__(self, source):
        """
        Source can be file path to json file or object loaded from the json file.
        """
        if isinstance(source, str):
            source = misc.load_json(source)
        self.W = source['W']
        self.PC = source['PC']
        self.QC = source['QC']
        self.groups = source['groups']
        # Convert coeffs arrays to numpy arrays.
        for group in self.groups:
            group['coeffs'] = np.array(group['coeffs']).reshape(-1, 6)

    def get_pos_vel(self, t: float):
        """
        Computes position and velocity.
        """
        t_pow = np.power(t, np.arange(6))
        t_pow_p = np.array([k * t_pow[k-1] if k > 0 else 0.0 for k in range(6)])

        v = np.zeros(3)
        vp = np.zeros(3)
        # Series terms have form c0 * t^\alpha * sin(c1 + c2*t + c3*t^2 + c4*t^3 + c5*t^4)
        for group in self.groups:
            coord = group['coord']
            alpha = group['alpha']
            coeffs = group['coeffs']

            a = coeffs[:,1].copy()
            ap = np.zeros(shape=a.shape)
            for k in range(1, 5):
                a += t_pow[k] * coeffs[:,k+1]
                ap += (k * t_pow[k-1]) * coeffs[:,k+1]
            c0_sin = np.sum(coeffs[:,0] * np.sin(a))
            c0_sin_p = np.sum(coeffs[:,0] * ap * np.cos(a))
            v[coord] += t_pow[alpha] * c0_sin
            vp[coord] += t_pow_p[alpha] * c0_sin + t_pow[alpha] * c0_sin_p

        # Compute spherical coordinates and their derivatives for the mean ecliptic of date.
        v[0] = v[0]*misc.ARCSEC + np.sum(self.W * t_pow[:5])
        v[1] *= misc.ARCSEC
        v[2] *= 0.9999999498265191
        vp[0] = vp[0]*misc.ARCSEC + np.sum(self.W * t_pow_p[:5])
        vp[1] *= misc.ARCSEC
        vp[2] *= 0.9999999498265191

        # Change to cartesian coordinates
        h = np.array([
            v[2]*np.cos(v[1])*np.cos(v[0]), 
            v[2]*np.cos(v[1])*np.sin(v[0]), 
            v[2]*np.sin(v[1])
        ])
        hp = np.array([
            (vp[2]*np.cos(v[1]) - vp[1]*v[2]*np.sin(v[1]))*np.cos(v[0]) - vp[0]*h[1],
            (vp[2]*np.cos(v[1]) - vp[1]*v[2]*np.sin(v[1]))*np.sin(v[0]) + vp[0]*h[0],
            vp[2]*np.sin(v[1]) + vp[1]*v[2]*np.cos(v[1])
        ])

        # Rotate from mean ecliptic and equinox of date to mean ecliptic and equinox of J2000.0.
        p = np.sum(self.PC * t_pow)
        q = np.sum(self.QC * t_pow)
        pp = np.sum(self.PC * t_pow_p)
        qp = np.sum(self.QC * t_pow_p)

        sc = np.sqrt(1.0 - p*p - q*q)
        pc1 = 1.0 - 2.0*p*p
        qc1 = 1.0 - 2.0*q*q
        pqp = pp*q + p*qp
        d2p = 2.0*p*pp + 2.0*q*qp
        pc2 = pp*sc - p*d2p/sc
        qc2 = qp*sc - q*d2p/sc

        pos = np.zeros(3)
        pos[0] = pc1*h[0] + 2.0*p*q*h[1] + 2.0*p*sc*h[2]
        pos[1] = 2.0*p*q*h[0] + qc1*h[1] - 2.0*q*sc*h[2]
        pos[2] = -2.0*p*sc*h[0] + 2.0*q*sc*h[1] + (pc1+qc1-1.0)*h[2]
        
        vel = np.zeros(3)
        vel[0] = pc1*hp[0] + 2.0*p*q*hp[1] + 2.0*p*sc*hp[2] 
        vel[0] += -4.0*p*pp*h[0] + 2.0*pqp*h[1] + 2.0*pc2*h[2]
        vel[1] = 2.0*p*q*hp[0] + qc1*hp[1] - 2.0*q*sc*hp[2] 
        vel[1] += 2.0*pqp*h[0] - 4.0*q*qp*h[1] - 2.0*qc2*h[2]
        vel[2] = -2.0*p*sc*hp[0] + 2.0*q*sc*hp[1] + (pc1+qc1-1.0)*hp[2] 
        vel[2] += -2.0*pc2*h[0] + 2.0*qc2*h[1] - 2.0*d2p*h[2]

        # Rotate from mean ecliptic and equinox of J2000.0 to mean equator and equinox of J2000.0.
        pos_equatorial = self.ECL_TO_EQU @ pos
        vel_equatorial = self.ECL_TO_EQU @ vel

        return pos_equatorial, vel_equatorial / 36525.0

def _test():
    # mpp02 = MPP02Ephemeris('./json/mpp02_llr_raw.json')
    mpp02 = MPP02Ephemeris('./json/mpp02_llr_truncated_medium.json')
    p, v = mpp02.get_pos_vel(-5.25)
    print(p.tolist(), v.tolist())

if __name__ == '__main__':
    _test()