import numpy as np

import tools

class VSOP87Ephemeris:
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
        # Covert from centuries to millenia:
        t = 0.1 * t
        pos = np.zeros(3)
        vel = np.zeros(3)
        for group in self.bodies[body_name]:
            coord_index = group['coord'] - 1
            alpha = group['alpha']
            coeffs = group['coeffs']

            sum_pos = np.sum(coeffs[:,0] * np.cos(coeffs[:,1] + t*coeffs[:,2]))
            sum_vel = -np.sum(coeffs[:,0] * coeffs[:,2] * np.sin(coeffs[:,1] + t*coeffs[:,2]))

            # Same without numpy vectorization:
            # sum_pos = 0.0
            # sum_vel = 0.0
            # for a, b, c in coeffs:
            #     sum_pos += a * np.cos(b + c*t)
            #     sum_vel -= a * c * np.sin(b + c*t)

            t_pow_alpha = np.power(t, alpha)
            t_pow_alpha_p = alpha * np.power(t, alpha-1) if alpha > 0 else 0.0
            pos[coord_index] += t_pow_alpha * sum_pos
            vel[coord_index] += t_pow_alpha_p * sum_pos + t_pow_alpha * sum_vel
        return self.matrix @ pos * self.AU, self.matrix @ vel * self.AU / 365250.0

if __name__ == '__main__':
    vsop87 = VSOP87Ephemeris('./json/vsop87a_truncated_7.json')
    p, v = vsop87.get_pos_vel('EARTH-MOON', 0.25)
    print(p, v)