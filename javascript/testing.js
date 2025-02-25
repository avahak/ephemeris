const fs = require('fs');

const VSOP87AEphemeris = require('./vsop87aEphemeris');
const MPP02Ephemeris = require('./mpp02Ephemeris');

function test() {
    let jsonVSOP, jsonMPP02;
    try {
        let data = fs.readFileSync('../json/vsop87a_raw.json', 'utf8');
        // let data = fs.readFileSync('../json/vsop87a_truncated_7.json', 'utf8');
        jsonVSOP = JSON.parse(data);

        data = fs.readFileSync('../json/mpp02_llr_raw.json', 'utf8');
        // data = fs.readFileSync('../json/mpp02_llr_truncated_7.json', 'utf8');
        jsonMPP02 = JSON.parse(data);
    } catch (err) {
        console.error('Error loading JSON file:', err);
        return;
    }

    // Test VSOP87AEphemeris
    const bodyNames = Object.keys(jsonVSOP['bodies']);
    let vsop87 = new VSOP87AEphemeris(jsonVSOP);
    let count = -0.5;
    let durationMs = 1000;
    let time0 = performance.now();
    while (performance.now() < time0 + durationMs) {
        let t = 20.0*Math.random() - 10.0;
        let bodyName = bodyNames[Math.floor(Math.random()*bodyNames.length)];
        const pv = vsop87.getPosVel(bodyName, t);
        count += 1;
        // console.log(count, bodyName, t, pv);
    }
    let pv = vsop87.getPosVel('EARTH-MOON', -5.25);
    console.log(`VSOP87A: getPosVel: ~${(count/durationMs*1000).toFixed(0)} calls/s.`);
    console.log("VSOP87A: getPosVel('EARTH-MOON', -5.25)", pv);


    // Test MPP02Ephemeris
    let mpp02 = new MPP02Ephemeris(jsonMPP02);
    count = -0.5;
    durationMs = 1000;
    time0 = performance.now();
    while (performance.now() < time0 + durationMs) {
        let t = 20.0*Math.random() - 10.0;
        const pv = mpp02.getPosVel(t);
        count += 1;
        // console.log(count, bodyName, t, pv);
    }
    pv = mpp02.getPosVel(-5.25);
    console.log(`MPP02: getPosVel: ~${(count/durationMs*1000).toFixed(0)} calls/s.`);
    console.log("MPP02: getPosVel(-5.25)", pv);
}

test();