// Copyright (c) 2017 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.

//! Algorithms for Lagrange interpolation.

use ::fields::Field;

pub struct LagrangeConstants<F: Field>(Vec<F::E>);

impl<F: Field> LagrangeConstants<F> 
{
    pub fn compute(point: &F::E, points: &[F::E], field: &F) -> LagrangeConstants<F> {
        let mut constants = Vec::with_capacity(points.len());
        for i in 0..points.len() {
            let xi = &points[i];
            let mut num = field.one();
            let mut denum = field.one();
            for j in 0..points.len() {
                if j != i {
                    let xj = &points[j];
                    num = field.mul(num, field.sub(xj, point));
                    denum = field.mul(denum, field.sub(xj, xi));
                }
            }
            let coef = field.mul(num, field.inv(denum));
            constants.push(coef);
        }
        
        LagrangeConstants(constants)
    }
    
    /// Note that care must be taken to provide the same `field` as the one used 
    /// for computing the constants!
    pub fn interpolate(&self, values: &[F::E], field: &F) -> F::E {
        let ref constants = self.0;
        assert_eq!(values.len(), constants.len());
        // compute weighted sum
        ::numtheory::weighted_sum(values, constants, field)
    }
}

/// Performs Lagrange interpolation at the specified point,
/// for a polynomial defined by `points` and `values`.
///
/// `points` and `values` are expected to be two arrays of the same size, containing
/// respectively the evaluation points (x) and the value of the polynomial at those point (p(x)).
pub fn lagrange_interpolation_at_point<F>(point: &F::E, points: &[F::E], values: &[F::E], field: &F) -> F::E
where F: Field, F::E: Clone
{
    assert_eq!(points.len(), values.len());
    let constants = LagrangeConstants::compute(point, points, field);
    constants.interpolate(values, field)
}

/// Performs Lagrange interpolation at the origin,
/// for a polynomial defined by `points` and `values`.
/// The result is the value of the polynomial at x=0.
///
/// `points` and `values` are expected to be two arrays of the same size, containing
/// respectively the evaluation points (x) and the value of the polynomial at those point (p(x)).
pub fn lagrange_interpolation_at_zero<F>(points: &[F::E], values: &[F::E], field: &F) -> F::E
where F: Field, F::E: Clone
{
    lagrange_interpolation_at_point(&field.zero(), points, values, field)
}

#[cfg(test)]
mod tests {
    
    use super::*;
    use ::fields::*;
    
    fn test_interpolation_from_constants<F>()
    where F: PrimeField + New<u32> + Encode<u32> + Decode<u32>, F::P: From<u32>, F::E: Clone
    {
        let ref field = F::new(17);

        let poly = field.encode_slice([4, 3, 2, 1]);
        let points = field.encode_slice([5, 6, 7, 8, 9]);
        
        let values = points.iter()
            .map(|point| ::numtheory::mod_evaluate_polynomial(&poly, point, field))
            .collect::<Vec<_>>();
        
        let constants = LagrangeConstants::compute(&field.zero(), &points, field);
        let value = constants.interpolate(&values, field);
        assert_eq!(field.decode(value), 4);
    }
    
    fn test_lagrange_interpolation_at_zero<F>() 
    where F: PrimeField + New<u32> + Encode<u32> + Decode<u32>, F::P: From<u32>, F::E: Clone
    {
        let ref field = F::new(17);

        let poly = field.encode_slice([4, 3, 2, 1]);
        let points = field.encode_slice([5, 6, 7, 8, 9]);
        
        let values = points.iter()
            .map(|point| ::numtheory::mod_evaluate_polynomial(&poly, point, field))
            .collect::<Vec<_>>();
            
        assert_eq!(field.decode_slice(&values), [7, 4, 7, 5, 4]);
        assert_eq!(field.decode(lagrange_interpolation_at_zero(&points, &values, field)), 4);
    }

    macro_rules! all_tests {
        ($field:ty) => {
            #[test] fn test_interpolation_from_constants() { super::test_interpolation_from_constants::<$field>(); }
            #[test] fn test_lagrange_interpolation_at_zero() { super::test_lagrange_interpolation_at_zero::<$field>(); }
        }
    }
    
    mod natural { all_tests!(::fields::NaturalPrimeField<i64>); }
    mod montgomery { all_tests!(::fields::MontgomeryField32); }
    #[cfg(feature="largefield_ramp")] mod large { all_tests!(::fields::LargePrimeField); }
    #[cfg(feature="largefield_gmp")] mod large { all_tests!(::fields::LargePrimeField); }
    
}
