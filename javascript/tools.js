/**
 * This can be replaced with math.js/equivalent calls.
 */
exports.vectorMult = (vector, c) => {
    const result = [...vector];
    for (let m = 0; m < vector.length; m++) 
        result[m] *= c;
    return result;
}

/**
 * This can be replaced with math.js/equivalent calls.
 */
exports.matrixMult = (matrix, vector) => {
    const result = new Array(vector.length).fill(0);
    for (let m = 0; m < matrix.length; m++) 
        for (let n = 0; n < matrix[m].length; n++)
            result[m] += matrix[m][n] * vector[n];
    return result;
}