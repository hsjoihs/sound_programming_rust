use sinc;
use std::f64::consts::PI;

#[allow(non_snake_case)]
pub fn safe_FIR_LPF(fe: f64, J: usize, b: &mut [f64], w: &mut [f64]) {
    assert_eq!(J % 2, 0);
    assert_eq!(J + 1, b.len());
    assert_eq!(J + 1, w.len());

    // k = m + J/2
    for (k, item) in b.iter_mut().enumerate() {
        *item = 2.0 * fe * sinc(2.0 * PI * fe * (k as f64 - (J / 2) as f64));
    }

    for (m, item) in b.iter_mut().enumerate() {
        *item *= w[m];
    }
}

#[allow(non_snake_case)]
pub fn safe_FIR_filtering(x: &[f64], y: &mut [f64], L: usize, b: &mut [f64], J: usize) {
    // check index here
    assert_eq!(J + 1, b.len());

    for n in 0..L {
        for m in 0..=J {
            if n >= m {
                y[n] += b[m] * x[n - m];
            }
        }
    }
}

// not tested
#[allow(non_snake_case)]
pub fn safe_FIR_HPF(fe: f64, J: usize, b: &mut [f64], w: &mut [f64]) {
    assert_eq!(J % 2, 0);
    assert_eq!(J + 1, b.len());
    assert_eq!(J + 1, w.len());
    let J = J as i32;
    let offset = J / 2;
    for m in -(J / 2)..J / 2 {
        b[(offset + m) as usize] = sinc(PI * m as f64) - 2.0 * fe * sinc(2.0 * PI * fe * m as f64);
    }

    for m in 0..(J + 1) as usize {
        b[m] *= w[m];
    }
}

// not tested
#[allow(non_snake_case)]
pub fn safe_FIR_BPF(fe1: f64, fe2: f64, J: usize, b: &mut [f64], w: &mut [f64]) {
    assert_eq!(J % 2, 0);
    assert_eq!(J + 1, b.len());
    assert_eq!(J + 1, w.len());
    let J = J as i32;
    let offset = J / 2;
    for m in -(J / 2)..J / 2 {
        b[(offset + m) as usize] = 2.0 * fe2 * sinc(2.0 * PI * fe2 * m as f64)
            - 2.0 * fe1 * sinc(2.0 * PI * fe1 * m as f64);
    }

    for m in 0..(J + 1) as usize {
        b[m] *= w[m];
    }
}

// not tested
#[allow(non_snake_case)]
pub fn safe_FIR_BEF(fe1: f64, fe2: f64, J: usize, b: &mut [f64], w: &mut [f64]) {
    assert_eq!(J % 2, 0);
    assert_eq!(J + 1, b.len());
    assert_eq!(J + 1, w.len());
    let J = J as i32;
    let offset = J / 2;
    for m in -(J / 2)..J / 2 {
        b[(offset + m) as usize] = sinc(PI * m as f64) - 2.0 * fe2 * sinc(2.0 * PI * fe2 * m as f64)
            + 2.0 * fe1 * sinc(2.0 * PI * fe1 * m as f64);
    }

    for m in 0..(J + 1) as usize {
        b[m] *= w[m];
    }
}

//not tested
#[allow(non_snake_case)]
pub fn safe_IIR_HPF(fc: f64, Q: f64, a: &mut [f64], b: &mut [f64]) {
    assert_eq!(3, a.len());
    assert_eq!(3, b.len());
    let fc = (PI * fc).tan() / (2.0 * PI);
    a[0] = 1.0 + 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc;
    a[1] = (8.0 * PI * PI * fc * fc - 2.0) / a[0];
    a[2] = (1.0 - 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc) / a[0];
    b[0] = 1.0 / a[0];
    b[1] = -2.0 / a[0];
    b[2] = 1.0 / a[0];

    a[0] = 1.0;
}

//not tested
#[allow(non_snake_case)]
pub fn safe_IIR_BPF(fc1: f64, fc2: f64, a: &mut [f64], b: &mut [f64]) {
    assert_eq!(3, a.len());
    assert_eq!(3, b.len());
    let fc1 = (PI * fc1).tan() / (2.0 * PI);
    let fc2 = (PI * fc2).tan() / (2.0 * PI);
    a[0] = 1.0 + 2.0 * PI * (fc2 - fc1) + 4.0 * PI * PI * fc1 * fc2;
    a[1] = (8.0 * PI * PI * fc1 * fc2 - 2.0) / a[0];
    a[2] = (1.0 - 2.0 * PI * (fc2 - fc1) + 4.0 * PI * PI * fc1 * fc2) / a[0];
    b[0] = 2.0 * PI * (fc2 - fc1) / a[0];
    b[1] = 0.0;
    b[2] = -2.0 * PI * (fc2 - fc1) / a[0];

    a[0] = 1.0;
}

//not tested
#[allow(non_snake_case)]
pub fn safe_IIR_BEF(fc1: f64, fc2: f64, a: &mut [f64], b: &mut [f64]) {
    assert_eq!(3, a.len());
    assert_eq!(3, b.len());
    let fc1 = (PI * fc1).tan() / (2.0 * PI);
    let fc2 = (PI * fc2).tan() / (2.0 * PI);

    a[0] = 1.0 + 2.0 * PI * (fc2 - fc1) + 4.0 * PI * PI * fc1 * fc2;
    a[1] = (8.0 * PI * PI * fc1 * fc2 - 2.0) / a[0];
    a[2] = (1.0 - 2.0 * PI * (fc2 - fc1) + 4.0 * PI * PI * fc1 * fc2) / a[0];
    b[0] = (4.0 * PI * PI * fc1 * fc2 + 1.0) / a[0];
    b[1] = (8.0 * PI * PI * fc1 * fc2 - 2.0) / a[0];
    b[2] = (4.0 * PI * PI * fc1 * fc2 + 1.0) / a[0];

    a[0] = 1.0;
}

//not tested
#[allow(non_snake_case)]
pub fn safe_IIR_notch(fc: f64, Q: f64, a: &mut [f64], b: &mut [f64]) {
    assert_eq!(3, a.len());
    assert_eq!(3, b.len());
    let fc = (PI * fc).tan() / (2.0 * PI);
    a[0] = 1.0 + 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc;
    a[1] = (8.0 * PI * PI * fc * fc - 2.0) / a[0];
    a[2] = (1.0 - 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc) / a[0];
    b[0] = (4.0 * PI * PI * fc * fc + 1.0) / a[0];
    b[1] = (8.0 * PI * PI * fc * fc - 2.0) / a[0];
    b[2] = (4.0 * PI * PI * fc * fc + 1.0) / a[0];

    a[0] = 1.0;
}

//not tested
#[allow(non_snake_case)]
pub fn safe_IIR_low_shelving(fc: f64, Q: f64, g: f64, a: &mut [f64], b: &mut [f64]) {
    assert_eq!(3, a.len());
    assert_eq!(3, b.len());
    let fc = (PI * fc).tan() / (2.0 * PI);
    a[0] = 1.0 + 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc;
    a[1] = (8.0 * PI * PI * fc * fc - 2.0) / a[0];
    a[2] = (1.0 - 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc) / a[0];
    b[0] =
        (1.0 + (1.0 + g).sqrt() * 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc * (1.0 + g)) / a[0];
    b[1] = (8.0 * PI * PI * fc * fc * (1.0 + g) - 2.0) / a[0];
    b[2] =
        (1.0 - (1.0 + g).sqrt() * 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc * (1.0 + g)) / a[0];

    a[0] = 1.0;
}

//not tested
#[allow(non_snake_case)]
pub fn safe_IIR_high_shelving(fc: f64, Q: f64, g: f64, a: &mut [f64], b: &mut [f64]) {
    assert_eq!(3, a.len());
    assert_eq!(3, b.len());
    let fc = (PI * fc).tan() / (2.0 * PI);
    a[0] = 1.0 + 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc;
    a[1] = (8.0 * PI * PI * fc * fc - 2.0) / a[0];
    a[2] = (1.0 - 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc) / a[0];
    b[0] = ((1.0 + g) + (1.0 + g).sqrt() * 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc) / a[0];
    b[1] = (8.0 * PI * PI * fc * fc - 2.0 * (1.0 + g)) / a[0];
    b[2] = ((1.0 + g) - (1.0 + g).sqrt() * 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc) / a[0];

    a[0] = 1.0;
}

//not tested
#[allow(non_snake_case)]
pub fn safe_IIR_peaking(fc: f64, Q: f64, g: f64, a: &mut [f64], b: &mut [f64]) {
    assert_eq!(3, a.len());
    assert_eq!(3, b.len());
    let fc = (PI * fc).tan() / (2.0 * PI);
    a[0] = 1.0 + 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc;
    a[1] = (8.0 * PI * PI * fc * fc - 2.0) / a[0];
    a[2] = (1.0 - 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc) / a[0];
    b[0] = (1.0 + 2.0 * PI * fc / Q * (1.0 + g) + 4.0 * PI * PI * fc * fc) / a[0];
    b[1] = (8.0 * PI * PI * fc * fc - 2.0) / a[0];
    b[2] = (1.0 - 2.0 * PI * fc / Q * (1.0 + g) + 4.0 * PI * PI * fc * fc) / a[0];

    a[0] = 1.0;
}

#[allow(non_snake_case)]
pub fn safe_IIR_LPF(fc: f64, Q: f64, a: &mut [f64], b: &mut [f64]) {
    assert_eq!(3, a.len());
    assert_eq!(3, b.len());
    let fc = (PI * fc).tan() / (2.0 * PI);

    a[0] = 1.0 + 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc;
    a[1] = (8.0 * PI * PI * fc * fc - 2.0) / a[0];
    a[2] = (1.0 - 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc) / a[0];
    b[0] = 4.0 * PI * PI * fc * fc / a[0];
    b[1] = 8.0 * PI * PI * fc * fc / a[0];
    b[2] = 4.0 * PI * PI * fc * fc / a[0];

    a[0] = 1.0;
}

#[allow(non_snake_case)]
pub fn safe_IIR_filtering(
    x: &[f64],
    y: &mut [f64],
    L: usize,
    a: &[f64],
    b: &[f64],
    I: usize,
    J: usize,
) {
    assert_eq!(J + 1, b.len());
    assert_eq!(I + 1, a.len());
    assert_eq!(L, x.len());
    assert_eq!(L, y.len());

    for n in 0..L {
        for m in 0..=J {
            if n >= m {
                y[n] += b[m] * x[n - m];
            }
        }
        for m in 1..=I {
            if n >= m {
                y[n] += -a[m] * y[n - m];
            }
        }
    }
}

#[allow(non_snake_case)]
pub fn safe_IIR_resonator(fc: f64, Q: f64, a: &mut [f64], b: &mut [f64]) {
    assert_eq!(3, a.len());
    assert_eq!(3, b.len());
    let fc = (PI * fc).tan() / (2.0 * PI);

    a[0] = 1.0 + 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc;
    a[1] = (8.0 * PI * PI * fc * fc - 2.0) / a[0];
    a[2] = (1.0 - 2.0 * PI * fc / Q + 4.0 * PI * PI * fc * fc) / a[0];
    b[0] = 2.0 * PI * fc / Q / a[0];
    b[1] = 0.0;
    b[2] = -2.0 * PI * fc / Q / a[0];

    a[0] = 1.0;
}
