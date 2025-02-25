import numpy as np

import tools

class VSOP87AEphemeris:
    """
    Computes planet positions and velocities w.r.t. the Sun using VSOP87A series.
    """
    AU = 149597870.691  # AU in km

    def __init__(self, source):
        """
        Here source can be file path to json file containing the series or
        the json object itself.
        """
        if isinstance(source, str):
            source = tools.load_json(source)
        self.bodies = source['bodies']
        # Reshape the coefficients from (n,) to (n/3,3)
        for body in self.bodies.values():
            for group in body:
                group['coeffs'] = np.array(group['coeffs']).reshape(-1, 3)
        self.matrix = np.array(source['matrix'])

    def get_pos_vel(self, body_name: str, t: float):
        """
        Returns position and velocity for given body and time in Julian centuries since J2000.0.
        Units are km, Julian day.
        """
        # Convert from centuries to millenia:
        t = 0.1 * t
        t_pow = np.power(t, np.arange(6))
        t_pow_p = np.array([k * t_pow[k-1] if k > 0 else 0.0 for k in range(6)])

        pos = np.zeros(3)
        vel = np.zeros(3)
        for group in self.bodies[body_name]:
            coord = group['coord']
            alpha = group['alpha']
            coeffs = group['coeffs']

            sum_pos = np.sum(coeffs[:,0] * np.cos(coeffs[:,1] + t*coeffs[:,2]))
            sum_vel = -np.sum(coeffs[:,0] * coeffs[:,2] * np.sin(coeffs[:,1] + t*coeffs[:,2]))

            pos[coord] += t_pow[alpha] * sum_pos
            vel[coord] += t_pow_p[alpha] * sum_pos + t_pow[alpha] * sum_vel
        return self.matrix @ pos * self.AU, self.matrix @ vel * self.AU / 365250.0

if __name__ == '__main__':
    vsop87 = VSOP87AEphemeris('./json/vsop87a_raw.json')
    # vsop87 = VSOP87AEphemeris('./json/vsop87a_truncated_7.json')
    p, v = vsop87.get_pos_vel('EARTH-MOON', -5.25)
    print(p.tolist(), v.tolist())