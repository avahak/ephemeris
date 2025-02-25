/**
 * 3x3 rotation matrix
 */
module.exports.rotationMatrix = (index, theta) => {
    const [c, s] = [Math.cos(theta), Math.sin(theta)];
    const [k1, k2] = [(index+1)%3, (index+2)%3];

    const rot = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    rot[k1][k1] = c;
    rot[k1][k2] = -s;
    rot[k2][k1] = s;
    rot[k2][k2] = c;
    return rot;
};

/**
 * This can be replaced with math.js/equivalent calls.
 */
module.exports.vectorMult = (vector, c) => {
    const result = [...vector];
    for (let m = 0; m < vector.length; m++) 
        result[m] *= c;
    return result;
};

/**
 * This can be replaced with math.js/equivalent calls.
 */
module.exports.matrixMult = (matrix, vector) => {
    const result = new Array(vector.length).fill(0);
    for (let m = 0; m < matrix.length; m++) 
        for (let n = 0; n < matrix[m].length; n++)
            result[m] += matrix[m][n] * vector[n];
    return result;
};