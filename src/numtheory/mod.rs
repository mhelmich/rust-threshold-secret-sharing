// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.

//! Various number theoretic utility functions used in the library.

pub mod fft;

use ::fields::Field;
use ::fields::Encode;

/// Euclidean GCD implementation (recursive). The first member of the returned
/// triplet is the GCD of `a` and `b`.
pub fn gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if b == 0 {
        (a, 1, 0)
    } else {
        let n = a / b;
        let c = a % b;
        let r = gcd(b, c);
        (r.0, r.2, r.1 - r.2 * n)
    }
}

#[test]
fn test_gcd() {
    assert_eq!(gcd(12, 16), (4, -1, 1));
}

pub fn binary_egcd(mut a: i64, mut b: i64) -> (i64, i64, i64) {
    
    // simple cases
    if a == 0 { return (b, 0, 1) }
    if b == 0 { return (a, 1, 0) }
    
    let mut u = 1;
    let mut v = 0;
    let mut s = 0;
    let mut t = 1;
    
    // find greatest power r of 2 dividing both a and b
    let mut r = 0;
    while (a | b) & 1 == 0 {
        a >>= 1;
        b >>= 1;
        r += 1;
    }
    
    let alpha = a;
    let beta = b;
    
    while a & 1 == 0 {
        a >>= 1;
        
        if (u | v) & 1 == 0 {
            u >>= 1;
            v >>= 1;
        } else {
            u = (u + beta) >> 1; 
            v = (v - alpha) >> 1;
        }
    }
    
    while a != b {
        if b & 1 == 0 {
            b >>= 1;
            if (s | t) & 1 == 0 {
                s >>= 1;
                t >>= 1;
            } else {
                s = (s + beta) >> 1;
                t = (t - alpha) >> 1;
            }
        } else if b < a {
            a = b;
            b = a;
            u = s;
            v = t;
            s = u;
            t = v;
        } else {
            b = b - a;
            s = s - u;
            t = t - v;
        }
    }
    
    (a << r, s, t)
}

#[test]
pub fn test_binary_egcd() {
    assert_eq!(binary_egcd(10, 4), (2, 1, -2));
}

/// Inverse of `k` in the *Zp* field defined by `prime`.
pub fn mod_inverse(k: i64, prime: i64) -> i64 {
    let k2 = k % prime;
    let r = if k2 < 0 {
        -gcd(prime, -k2).2
    } else {
        gcd(prime, k2).2
    };
    (prime + r) % prime
}
// pub fn mod_inverse(k: i64, prime: i64) -> i64 {
//     let k2 = k % prime;
//     let r = if k2 < 0 {
//         -binary_egcd(prime, -k2).2
//     } else {
//         binary_egcd(prime, k2).2
//     };
//     (prime + r) % prime
// }

#[test]
fn test_mod_inverse() {
    assert_eq!(mod_inverse(3, 7), 5);
}


/// `x` to the power of `e` in the *Zp* field defined by `prime`.
pub fn mod_pow(mut x: i64, mut e: u32, prime: i64) -> i64 {
    let mut acc = 1;
    while e > 0 {
        if e % 2 == 0 {
            // even
            // no-op
        } else {
            // odd
            acc = (acc * x) % prime;
        }
        x = (x * x) % prime; // waste one of these by having it here but code is simpler (tiny bit)
        e = e >> 1;
    }
    acc
}

#[test]
fn test_mod_pow() {
    assert_eq!(mod_pow(2, 0, 17), 1);
    assert_eq!(mod_pow(2, 3, 17), 8);
    assert_eq!(mod_pow(2, 6, 17), 13);

    assert_eq!(mod_pow(-3, 0, 17), 1);
    assert_eq!(mod_pow(-3, 1, 17), -3);
    assert_eq!(mod_pow(-3, 15, 17), -6);
}

/// Performs a Lagrange interpolation in field Zp at the origin
/// for a polynomial defined by `points` and `values`.
///
/// `points` and `values` are expected to be two arrays of the same size, containing
/// respectively the evaluation points (x) and the value of the polynomial at those point (p(x)).
///
/// The result is the value of the polynomial at x=0. It is also its zero-degree coefficient.
///
/// This is obviously less general than `newton_interpolation_general` as we
/// only get a single value, but it is much faster.
pub fn lagrange_interpolation_at_zero<F>(points: &[F::E], values: &[F::E], field: &F) -> F::E
where F: Field, F::E: Copy
{
    // TODO coef computations could be reused
    
    assert_eq!(points.len(), values.len());
    // Lagrange interpolation for point 0
    let mut acc = field.zero();
    for i in 0..values.len() {
        // compute Lagrange coefficient
        let xi = points[i];
        let mut num = field.one();
        let mut denum = field.one();
        for j in 0..values.len() {
            if j != i {
                let xj = points[j];
                num = field.mul(num, xj);
                denum = field.mul(denum, field.sub(xj, xi));
            }
        }
        let coef = field.mul(num, field.inv(denum));
        // update sum
        let yi = values[i];
        acc = field.add(acc, field.mul(yi, coef));
    }
    acc
}

/// Holds together points and Newton-interpolated coefficients for fast evaluation.
pub struct NewtonPolynomial<F> 
where F: Field
{
    points: Vec<F::E>,
    coefficients: Vec<F::E>,
}


/// General case for Newton interpolation in field Zp.
///
/// Given enough `points` (x) and `values` (p(x)), find the coefficients for `p`.
pub fn newton_interpolation_general<F>(points: &[F::E],
                                       values: &[F::E],
                                       field: &F)
                                    -> NewtonPolynomial<F> 
where F: Field, F::E: Copy
{
    let coefficients = compute_newton_coefficients(points, values, field);
    NewtonPolynomial {
        points: points.to_vec(),
        coefficients: coefficients,
    }
}

#[test]
fn test_newton_interpolation_general() {
    let prime = 17;
    let ref field = ::fields::natural::NaturalPrimeField(prime);

    let poly = [1, 2, 3, 4];
    let points: Vec<i64> = vec![5, 6, 7, 8, 9];
    let values: Vec<i64> = points.iter()
        .map(|&point| mod_evaluate_polynomial(&poly, point, field))
        .collect();
    assert_eq!(values, vec![8, 16, 4, 13, 16]);

    let recovered_poly = newton_interpolation_general(&points, &values, field);
    let recovered_values: Vec<i64> = points.iter()
        .map(|&point| newton_evaluate(&recovered_poly, point, field))
        .collect();
    assert_eq!(recovered_values, values);

    assert_eq!(newton_evaluate(&recovered_poly, 10, field), 3);
    assert_eq!(newton_evaluate(&recovered_poly, 11, field), -2);
    assert_eq!(newton_evaluate(&recovered_poly, 12, field), 8);
}

pub fn newton_evaluate<F>(poly: &NewtonPolynomial<F>, point: F::E, field: &F) -> F::E 
where F: Field, F::E: Copy
{
    // compute Newton points
    let mut newton_points = vec![field.one()];
    for i in 0..poly.points.len() - 1 {
        let diff = field.sub(point, poly.points[i]);
        let product = field.mul(newton_points[i], diff);
        newton_points.push(product);
    }
    let ref newton_coefs = poly.coefficients;
    // sum up
    newton_coefs.iter()
        .zip(newton_points)
        .map(|(&coef, point)| field.mul(coef, point))
        .fold(field.zero(), |a, b| field.add(a, b))
}

fn compute_newton_coefficients<F>(points: &[F::E], values: &[F::E], field: &F) -> Vec<F::E> 
where F: Field, F::E: Copy
{
    assert_eq!(points.len(), values.len());

    let mut store: Vec<(usize, usize, F::E)> = 
        values.iter().enumerate()
        .map(|(index, &value)| (index, index, value))
        .collect();

    for j in 1..store.len() {
        for i in (j..store.len()).rev() {
            let index_lower = store[i - 1].0;
            let index_upper = store[i].1;

            let point_lower = points[index_lower];
            let point_upper = points[index_upper];
            let point_diff = field.sub(point_upper, point_lower);
            let point_diff_inverse = field.inv(point_diff);

            let coef_lower = store[i - 1].2;
            let coef_upper = store[i].2;
            let coef_diff = field.sub(coef_upper, coef_lower);

            let fraction = field.mul(coef_diff, point_diff_inverse);

            store[i] = (index_lower, index_upper, fraction);
        }
    }

    store.iter().map(|&(_, _, v)| v).collect()
}

#[test]
fn test_compute_newton_coefficients() {
    let points = vec![5, 6, 7, 8, 9];
    let values = vec![8, 16, 4, 13, 16];
    let field = ::fields::natural::NaturalPrimeField(17);

    let coefficients = compute_newton_coefficients(&points, &values, &field);
    assert_eq!(coefficients, vec![8, 8, -10, 4, 0]);
}

/// Map `values` from `[-n/2, n/2)` to `[0, n)`.
pub fn positivise(values: &[i64], n: i64) -> Vec<i64> {
    values.iter()
        .map(|&value| if value < 0 { value + n } else { value })
        .collect()
}

/// Evaluate polynomial given by `coefficients` at `point`.
///
/// Current implementation uses Horner's method.
pub fn mod_evaluate_polynomial<F>(coefficients: &[F::E], point: F::E, field: &F) -> F::E
where F: Field, F: Encode<u32>, F::E: Copy
{
    // evaluate using Horner's rule
    //  - to combine with fold we consider the coefficients in reverse order
    coefficients.iter().rev()
        .fold(field.zero(), |partial, &coef| field.add(field.mul(partial, point), coef))
}

#[test]
fn test_mod_evaluate_polynomial() {
    let poly = vec![1, 2, 3, 4, 5, 6];
    let point = 5;
    let field = ::fields::natural::NaturalPrimeField(17);
    assert_eq!(mod_evaluate_polynomial(&poly, point, &field), 4);
}
