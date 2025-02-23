"""
Calculates position and velocity of the Moon w.r.t. the Earth using ELP/MPP02.
Can use LLR or DE405 fit data and their truncated versions.

ELP/MMP02 article: Lunar solution ELP, version ELP/MPP02, Jean Chapront and Gerard Francou. 
"""
import numpy as np

import tools

class MPP02Ephemeris:
    """
    Computes position and velocity of the Moon using ELP/MPP02.
    Units: km for position, km/Julian day for velocity.
    """
    PC = np.array([0.10180391e-4, 0.47020439e-6, -0.5417367e-9, -0.2507948e-11, 0.463486e-14])
    QC = np.array([-0.113469002e-3, 0.12372674e-6, 0.1265417e-8, -0.1371808e-11, -0.320334e-14])
    EP = 0.40909280422232897    # ~23.44 deg

    def __init__(self, source):
        """
        Source can be file path to json file or object loaded from the json file.
        """
        if isinstance(source, str):
            source = tools.load_json(source)
        self.W = source['W']
        self.groups = source['groups']
        # Convert coeffs arrays to numpy arrays.
        for group in self.groups:
            group['coeffs'] = np.array(group['coeffs']).reshape(-1, 6)

    def get_pos(self, t: float):
        """
        Computes position only.
        """
        t_pow = np.power(t, np.arange(5))

        v = np.zeros(3)
        # Series terms have form c0 * t^\alpha * sin(c1 + c2*t + c3*t^2 + c4*t^3 + c5*t^4)
        for group in self.groups:
            coord = group['coord']
            alpha = group['alpha']
            coeffs = group['coeffs']

            a = coeffs[:,1].copy()
            for k in range(1, 5):
                a += coeffs[:,k+1] * t_pow[k]
            v[coord] += t_pow[alpha] * np.sum(coeffs[:,0] * np.sin(a))

            # Same without numpy vectorization:
            # for c in coeffs:
            #     a = c[1]
            #     for k in range(1, 5):
            #         a += c[k+1] * t_pow[k]
            #     v[coord] += c[0] * t_pow[alpha] * np.sin(a)

        # Compute spherical coordinates and their derivatives for the mean ecliptic of date.
        v[0] = v[0] * tools.ARCSEC + np.sum(self.W * t_pow)
        v[1] = v[1] * tools.ARCSEC
        v[2] = v[2] * 0.9999999498265191

        # Change to cartesian coordinates
        h = np.array([
            v[2]*np.cos(v[1])*np.cos(v[0]), 
            v[2]*np.cos(v[1])*np.sin(v[0]), 
            v[2]*np.sin(v[1])
        ])

        # Next, rotate from mean ecliptic and equinox of date to mean ecliptic and equinox of J2000.0.
        p = t * sum(self.PC * t_pow)
        q = t * sum(self.QC * t_pow)

        sc = np.sqrt(max(1.0 - p*p - q*q, 0.0))
        pc1 = 1.0 - 2.0*p*p
        qc1 = 1.0 - 2.0*q*q

        pos = np.zeros(3)
        pos[0] = pc1*h[0] + 2.0*p*q*h[1] + 2.0*p*sc*h[2]
        pos[1] = 2.0*p*q*h[0] + qc1*h[1] - 2.0*q*sc*h[2]
        pos[2] = -2.0*p*sc*h[0] + 2.0*q*sc*h[1] + (pc1+qc1-1.0)*h[2]

        # Finally, rotate from mean ecliptic of J2000.0 to mean equator and equinox of J2000.0.
        pos_equatorial = np.array([
            pos[0], 
            pos[1]*np.cos(self.EP) - pos[2]*np.sin(self.EP), 
            pos[1]*np.sin(self.EP) + pos[2]*np.cos(self.EP)
        ])

        return pos_equatorial

    def get_pos_vel(self, t: float):
        """
        Computes position and velocity.
        """
        t_pow = np.power(t, np.arange(5))

        v = np.zeros(3)
        vp = np.zeros(3)
        # Series terms have form c0 * t^\alpha * sin(c1 + c2*t + c3*t^2 + c4*t^3 + c5*t^4)
        for group in self.groups:
            coord = group['coord']
            alpha = group['alpha']
            coeffs = group['coeffs']
            t_pow_alpha_p = alpha * t_pow[alpha-1] if alpha > 0 else 0.0

            a = coeffs[:,1].copy()
            ap = np.zeros(shape=a.shape)
            for k in range(1, 5):
                a += t_pow[k] * coeffs[:,k+1]
                ap += (k * t_pow[k-1]) * coeffs[:,k+1]
            c0_sin = np.sum(coeffs[:,0] * np.sin(a))
            c0_sin_p = np.sum(coeffs[:,0] * ap * np.cos(a))
            v[coord] += t_pow[alpha] * c0_sin
            vp[coord] += t_pow_alpha_p * c0_sin + t_pow[alpha] * c0_sin_p

            # Same without numpy vectorization:
            # for c in coeffs:
            #     a = c[1]
            #     ap = 0.0
            #     for k in range(1, 5):
            #         a += c[k+1] * t_pow[k]
            #         ap += k * c[k+1] * t_pow[k-1]
            #     v[coord] += c[0] * t_pow[alpha] * np.sin(a)
            #     vp[coord] += c[0] * (t_pow_alpha_p * np.sin(a)  +  t_pow[alpha] * ap * np.cos(a))

        # Compute spherical coordinates and their derivatives for the mean ecliptic of date.
        v[0] = v[0] * tools.ARCSEC + np.sum(self.W * t_pow)
        v[1] = v[1] * tools.ARCSEC
        v[2] = v[2] * 0.9999999498265191
        vp[0] = vp[0] * tools.ARCSEC + sum([k*self.W[k]*t_pow[k-1] for k in range(1, 5)])
        vp[1] = vp[1] * tools.ARCSEC

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

        # Next, rotate from mean ecliptic and equinox of date to mean ecliptic and equinox of J2000.0.
        p = t * sum(self.PC * t_pow)
        q = t * sum(self.QC * t_pow)
        pp = sum([(k+1.0)*t_pow[k]*self.PC[k] for k in range(5)])
        qp = sum([(k+1.0)*t_pow[k]*self.QC[k] for k in range(5)])

        sc = np.sqrt(max(1.0 - p*p - q*q, 0.0))
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

        # Finally, rotate from mean ecliptic of J2000.0 to mean equator and equinox of J2000.0.
        pos_equatorial = np.array([
            pos[0], 
            pos[1]*np.cos(self.EP) - pos[2]*np.sin(self.EP), 
            pos[1]*np.sin(self.EP) + pos[2]*np.cos(self.EP)
        ])
        vel_equatorial = np.array([
            vel[0], 
            vel[1]*np.cos(self.EP) - vel[2]*np.sin(self.EP), 
            vel[1]*np.sin(self.EP) + vel[2]*np.cos(self.EP)
        ])

        return pos_equatorial, vel_equatorial / 36525.0

# if __name__ == '__main__':
#     pass