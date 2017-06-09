// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.

//! Algorithms for Newton interpolation.

use std::borrow::Borrow;
use ::fields::Field;

/// Holds together points and Newton-interpolated coefficients for fast evaluation.
pub struct NewtonPolynomial<F> 
where F: Field
{
    points: Vec<F::E>,
    coefficients: Vec<F::E>,
}

impl<F> NewtonPolynomial<F>
where F: Field, F::E: Clone
{
    /// General case for Newton interpolation in field Zp.
    ///
    /// Given enough `points` (x) and `values` (p(x)), find the coefficients for `p`.
    pub fn compute(points: &[F::E], values: &[F::E], field: &F) -> NewtonPolynomial<F>
    {
        assert_eq!(points.len(), values.len());

        // compute Newton coefficients
        
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

        let coefficients = store.into_iter().map(|(_, _, v)| v).collect();
        
        NewtonPolynomial {
            points: points.to_vec(),
            coefficients: coefficients,
        }
    }
    
    /// Note that care must be taken to ensure that the `field` is used here is the same 
    /// as the one used in `compute`.
    pub fn evaluate<P: Borrow<F::E>>(&self, point: P, field: &F) -> F::E 
    {
        // compute Newton points
        let mut newton_points = Vec::with_capacity(self.points.len() - 1);
        newton_points.push(field.one());
        for i in 0..self.points.len() - 1 {
            let diff = field.sub(point.borrow(), &self.points[i]);
            let product = field.mul(&newton_points[i], diff);
            newton_points.push(product);
        }
        
        ::numtheory::weighted_sum(&self.coefficients, &newton_points, field)
    }
}

pub fn newton_interpolation_at_point<F>(point: &F::E, points: &[F::E], values: &[F::E], field: &F) -> F::E
where F: Field, F::E: Clone
{
    assert_eq!(points.len(), values.len());
    let poly = NewtonPolynomial::compute(points, values, field);
    poly.evaluate(point, field)
}


#[cfg(test)]
mod tests {

    use super::*;
    use ::numtheory;
    use ::fields::*;

    fn test_newton_interpolation_general<F>() 
    where F: PrimeField + Encode<u32> + Decode<u32>, F::P: From<u32>, F::E: Clone
    {
        let ref field = F::new(17.into());

        let poly = field.encode_slice([1, 2, 3, 4]);
        let points = field.encode_slice([5, 6, 7, 8, 9]);
            
        let values = points.iter()
            .map(|point| numtheory::mod_evaluate_polynomial(&poly, point, field))
            .collect::<Vec<_>>();
        assert_eq!(field.decode_slice(&values), vec![8, 16, 4, 13, 16]);

        let recovered_poly = NewtonPolynomial::compute(&points, &values, field);
        let recovered_values = points.iter()
            .map(|point| recovered_poly.evaluate(point, field))
            .collect::<Vec<_>>();
        assert_eq!(field.decode_slice(recovered_values), field.decode_slice(values));

        assert_eq!(field.decode(recovered_poly.evaluate(field.encode(10), field)), 3);
        assert_eq!(field.decode(recovered_poly.evaluate(field.encode(11), field)), 15);
        assert_eq!(field.decode(recovered_poly.evaluate(field.encode(12), field)), 8);
    }

    fn test_compute_newton_coefficients<F>()
    where F: PrimeField + Encode<u32> + Decode<u32>, F::P: From<u32>, F::E: Clone
    {
        let ref field = F::new(17.into());
        
        let points = field.encode_slice([5, 6, 7, 8, 9]);
        let values = field.encode_slice([8, 16, 4, 13, 16]);

        let poly = NewtonPolynomial::compute(&points, &values, field);
        assert_eq!(field.decode_slice(poly.coefficients), vec![8, 8, 7, 4, 0]);
    }
    
    macro_rules! all_tests {
        ($field:ty) => {
            use super::*;
            #[test] fn test_newton_interpolation_general() { super::test_newton_interpolation_general::<$field>(); }
            #[test] fn test_compute_newton_coefficients() { super::test_compute_newton_coefficients::<$field>(); }
        } 
    }
    
    mod natural { all_tests!(::fields::NaturalPrimeField<i64>); }
    mod montgomery { all_tests!(::fields::MontgomeryField32); }
    #[cfg(feature="largefield")] mod large { all_tests!(::fields::LargePrimeField); }
    
}
