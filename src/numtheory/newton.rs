// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.

//! TODO

use ::fields::Field;

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
where F: Field, F::E: Clone
{
    let coefficients = compute_newton_coefficients(points, values, field);
    NewtonPolynomial {
        points: points.to_vec(),
        coefficients: coefficients,
    }
}

pub fn newton_evaluate<F>(poly: &NewtonPolynomial<F>, point: F::E, field: &F) -> F::E 
where F: Field, F::E: Clone
{
    // compute Newton points
    let mut newton_points = vec![field.one()];
    for i in 0..poly.points.len() - 1 {
        let diff = field.sub(&point, &poly.points[i]);
        let product = field.mul(&newton_points[i], diff);
        newton_points.push(product);
    }
    let ref newton_coefs = poly.coefficients;
    // sum up
    newton_coefs.iter()
        .zip(newton_points)
        .map(|(coef, point)| field.mul(coef, point))
        .fold(field.zero(), |a, b| field.add(a, b))
}

fn compute_newton_coefficients<F>(points: &[F::E], values: &[F::E], field: &F) -> Vec<F::E> 
where F: Field, F::E: Clone
{
    assert_eq!(points.len(), values.len());

    let mut store: Vec<(usize, usize, F::E)> = 
        values.iter().enumerate()
        .map(|(index, value)| (index, index, value.clone()))
        .collect();

    for j in 1..store.len() {
        for i in (j..store.len()).rev() {
            let index_lower = store[i - 1].0;
            let index_upper = store[i].1;

            let point_lower = &points[index_lower];
            let point_upper = &points[index_upper];
            let point_diff = field.sub(point_upper, point_lower);
            let point_diff_inverse = field.inv(point_diff);

            let coef_lower = store[i - 1].2.clone();
            let coef_upper = store[i].2.clone();
            let coef_diff = field.sub(coef_upper, coef_lower);

            let fraction = field.mul(coef_diff, point_diff_inverse);

            store[i] = (index_lower, index_upper, fraction);
        }
    }

    store.into_iter().map(|(_, _, v)| v).collect()
}


#[cfg(test)]
mod tests {

    use super::*;
    use numtheory;
    use fields;

    #[test]
    fn test_newton_interpolation_general() {
        let prime = 17;
        let ref field = fields::NaturalPrimeField(prime);

        let poly = [1, 2, 3, 4];
        let points: Vec<i64> = vec![5, 6, 7, 8, 9];
        let values: Vec<i64> = points.iter()
            .map(|&point| numtheory::mod_evaluate_polynomial(&poly, point, field))
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

    #[test]
    fn test_compute_newton_coefficients() {
        let points = vec![5, 6, 7, 8, 9];
        let values = vec![8, 16, 4, 13, 16];
        let field = fields::NaturalPrimeField(17);

        let coefficients = compute_newton_coefficients(&points, &values, &field);
        assert_eq!(coefficients, vec![8, 8, -10, 4, 0]);
    }
    
}
