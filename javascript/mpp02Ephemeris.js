const { rotationMatrix, matrixMult, vectorMult } = require('./tools');

/**
 * Computes position and velocity of the Moon using ELP/MPP02.
 * Units: km for position, km/Julian day for velocity.
 */
class MPP02Ephemeris {
    PC = new Float64Array([0, 0.10180391e-4, 0.47020439e-6, -0.5417367e-9, -0.2507948e-11, 0.463486e-14]);
    QC = new Float64Array([0, -0.113469002e-3, 0.12372674e-6, 0.1265417e-8, -0.1371808e-11, -0.320334e-14]);
    ARCSEC = Math.PI / 180 / 3600;
    ECL_TO_EQU = rotationMatrix(0, 84381.448 * this.ARCSEC);    // angle ~23.44 deg

    /**
     * Here source has to be object containing the (possibly truncated) ELP/MPP02 series
     * loaded from one of the json files.
     */
    constructor(source) {
        this.W = source['W'];
        this.groups = source['groups'];
    }

    /**
     * Returns position and velocity for the Moon given time in Julian centuries since J2000.0. 
     * Units are km, Julian day.
     */
    getPosVel(t) {
        const tPow = Array.from({ length: 6 }, (_, k) => Math.pow(t, k));
        const tPowP = Array.from({ length: 6 }, (_, k) => k > 0 ? k * tPow[k-1] : 0);

        let v = [0, 0, 0];
        let vp = [0, 0, 0];
        // Series terms have form c0 * t^\alpha * sin(c1 + c2*t + c3*t^2 + c4*t^3 + c5*t^4)
        this.groups.forEach((group) => {
            const coord = group['coord'];
            const alpha = group['alpha'];
            const coeffs = group['coeffs'];

            for (let k6 = 0; k6 < coeffs.length; k6 += 6) {
                let a = 0;
                let ap = 0;
                for (let j = 0; j < 5; j++) {
                    const c = coeffs[k6 + j + 1];
                    a += c * tPow[j];
                    ap += c * tPowP[j];
                }
                v[coord] += coeffs[k6 + 0] * tPow[alpha] * Math.sin(a);
                vp[coord] += coeffs[k6 + 0] * (tPowP[alpha] * Math.sin(a) + tPow[alpha] * ap * Math.cos(a));
            }
        });

        // Compute spherical coordinates and their derivatives for the mean ecliptic of date.
        v[0] = v[0]*this.ARCSEC + this.W.reduce((acc, value, k) => acc + value * tPow[k], 0);
        v[1] *= this.ARCSEC;
        v[2] *= 0.9999999498265191;
        vp[0] = vp[0]*this.ARCSEC + this.W.reduce((acc, value, k) => acc + value * tPowP[k], 0);
        vp[1] *= this.ARCSEC;
        vp[2] *= 0.9999999498265191;

        // Change to cartesian coordinates
        const h = [
            v[2]*Math.cos(v[1])*Math.cos(v[0]), 
            v[2]*Math.cos(v[1])*Math.sin(v[0]), 
            v[2]*Math.sin(v[1])
        ];
        const hp = [
            (vp[2]*Math.cos(v[1]) - vp[1]*v[2]*Math.sin(v[1]))*Math.cos(v[0]) - vp[0]*h[1],
            (vp[2]*Math.cos(v[1]) - vp[1]*v[2]*Math.sin(v[1]))*Math.sin(v[0]) + vp[0]*h[0],
            vp[2]*Math.sin(v[1]) + vp[1]*v[2]*Math.cos(v[1])
        ];

        // Next, rotate from mean ecliptic and equinox of date to mean ecliptic and equinox of J2000.0.
        const p = this.PC.reduce((acc, value, k) => acc + value * tPow[k], 0);
        const q = this.QC.reduce((acc, value, k) => acc + value * tPow[k], 0);
        const pp = this.PC.reduce((acc, value, k) => acc + value * tPowP[k], 0);
        const qp = this.QC.reduce((acc, value, k) => acc + value * tPowP[k], 0);

        const sc = Math.sqrt(1 - p*p - q*q);
        const pc1 = 1 - 2*p*p;
        const qc1 = 1 - 2*q*q;
        const pqp = pp*q + p*qp;
        const d2p = 2*p*pp + 2*q*qp;
        const pc2 = pp*sc - p*d2p/sc;
        const qc2 = qp*sc - q*d2p/sc;

        const pos = [0, 0, 0];
        pos[0] = pc1*h[0] + 2*p*q*h[1] + 2*p*sc*h[2];
        pos[1] = 2*p*q*h[0] + qc1*h[1] - 2*q*sc*h[2];
        pos[2] = -2*p*sc*h[0] + 2*q*sc*h[1] + (pc1+qc1-1)*h[2];
        
        const vel = [0, 0, 0];
        vel[0] = pc1*hp[0] + 2*p*q*hp[1] + 2*p*sc*hp[2];
        vel[0] += -4*p*pp*h[0] + 2*pqp*h[1] + 2*pc2*h[2];
        vel[1] = 2*p*q*hp[0] + qc1*hp[1] - 2*q*sc*hp[2];
        vel[1] += 2*pqp*h[0] - 4*q*qp*h[1] - 2*qc2*h[2];
        vel[2] = -2*p*sc*hp[0] + 2*q*sc*hp[1] + (pc1+qc1-1)*hp[2];
        vel[2] += -2*pc2*h[0] + 2*qc2*h[1] - 2*d2p*h[2];

        // Finally, rotate from mean ecliptic and equinox of J2000.0 to mean equator and equinox of J2000.0.
        const posEquatorial = matrixMult(this.ECL_TO_EQU, pos);
        const velEquatorial = matrixMult(this.ECL_TO_EQU, vel);

        return { 'pos': posEquatorial, 'vel': vectorMult(velEquatorial, 1/36525) };
    }
}

module.exports = MPP02Ephemeris;