/**
 * Sample code loading json files and computing positions and velocities
 * with VSOP87AEphemeris and MPP02Ephemeris. Execution time is estimated also.
 */

const fs = require('fs');

const VSOP87AEphemeris = require('./vsop87aEphemeris');
const MPP02Ephemeris = require('./mpp02Ephemeris');

function test(size) {
    // Load JSON files
    let jsonVSOP, jsonMPP02;
    try {
        // let data = fs.readFileSync('../json/vsop87a_raw.json', 'utf8');
        let data = fs.readFileSync(`../json/vsop87a_truncated_${size}.json`, 'utf8');
        jsonVSOP = JSON.parse(data);

        // data = fs.readFileSync('../json/mpp02_llr_raw.json', 'utf8');
        data = fs.readFileSync(`../json/mpp02_llr_truncated_${size}.json`, 'utf8');
        jsonMPP02 = JSON.parse(data);
    } catch (err) {
        console.error('Error loading JSON file:', err);
        return;
    }

    // Test VSOP87AEphemeris
    const bodyNames = Object.keys(jsonVSOP['bodies']);
    let vsop87 = new VSOP87AEphemeris(jsonVSOP);
    let count = 0;
    let durationMs = 1500;
    let time0 = performance.now();
    let time = time0;
    while (time < time0 + durationMs) {
        let t = 60.0*Math.random() - 30.0;
        let bodyName = bodyNames[Math.floor(Math.random()*bodyNames.length)];
        const pv = vsop87.getPosVel(bodyName, t);
        time = performance.now();
        count += 1;
        // console.log(count, bodyName, t, pv);
    }
    let pv = vsop87.getPosVel('EARTH-MOON', -5.25);
    console.log(`VSOP87A: getPosVel: ~${(count/(time-time0)*1000).toFixed(0)} calls/s.`);
    console.log("VSOP87A: getPosVel('EARTH-MOON', -5.25)", pv);


    // Test MPP02Ephemeris
    let mpp02 = new MPP02Ephemeris(jsonMPP02);
    count = 0;
    durationMs = 1500;
    time0 = performance.now();
    time = time0;
    while (time < time0 + durationMs) {
        let t = 60.0*Math.random() - 30.0;
        const pv = mpp02.getPosVel(t);
        time = performance.now();
        count += 1;
        // console.log(count, t, pv);
    }
    pv = mpp02.getPosVel(-5.25);
    console.log(`MPP02: getPosVel: ~${(count/(time-time0)*1000).toFixed(0)} calls/s.`);
    console.log("MPP02: getPosVel(-5.25)", pv);
}

// test('small');
test('medium');
// test('large');