// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.

//! TODO

use ::fields::Field;

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
where F: Field, F::E: Clone
{
    // TODO coef computations could be reused
    
    assert_eq!(points.len(), values.len());
    // Lagrange interpolation for point 0
    let mut acc = field.zero();
    for i in 0..values.len() {
        // compute Lagrange coefficient
        let xi = &points[i];
        let mut num = field.one();
        let mut denum = field.one();
        for j in 0..values.len() {
            if j != i {
                let xj = &points[j];
                num = field.mul(num, xj);
                denum = field.mul(denum, field.sub(xj, xi));
            }
        }
        let coef = field.mul(num, field.inv(denum));
        // update sum
        let yi = &values[i];
        acc = field.add(acc, field.mul(yi, coef));
    }
    acc
}

// TODO
// #[test]
// fn lagrange_interpolation_at_zero() {
//     let prime = 17;
//     let ref field = ::fields::natural::NaturalPrimeField(prime);
// 
//     let poly = [1, 2, 3, 4];
//     let points: Vec<i64> = vec![5, 6, 7, 8, 9];
//     let values: Vec<i64> = points.iter()
//         .map(|&point| ::numtheory::mod_evaluate_polynomial(&poly, point, field))
//         .collect();
//     assert_eq!(values, vec![8, 16, 4, 13, 16]);
// 
//     let recovered_poly = lagrange_interpolation_at_zero(&points, &values, field);
//     let recovered_values: Vec<i64> = points.iter()
//         .map(|&point| newton_evaluate(&recovered_poly, point, field))
//         .collect();
//     assert_eq!(recovered_values, values);
// 
//     assert_eq!(newton_evaluate(&recovered_poly, 10, field), 3);
//     assert_eq!(newton_evaluate(&recovered_poly, 11, field), -2);
//     assert_eq!(newton_evaluate(&recovered_poly, 12, field), 8);
// }
