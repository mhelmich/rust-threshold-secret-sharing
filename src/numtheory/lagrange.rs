// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.

//! Algorithms for Lagrange interpolation.

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

#[cfg(test)]
mod tests {
    
    use super::*;
    use ::fields::*;
    
    fn test_lagrange_interpolation_at_zero<F>() 
    where F: PrimeField + Encode<u32> + Decode<u32>, F::P: From<u32>, F::E: Clone
    {
        let ref field = F::new(17.into());

        let poly = field.encode_slice([4, 3, 2, 1]);
        let points = field.encode_slice([5, 6, 7, 8, 9]);
        
        let values = points.iter()
            .map(|point| ::numtheory::mod_evaluate_polynomial(&poly, point.clone(), field)) // TODO no need to clone
            .collect::<Vec<_>>();
            
        assert_eq!(field.decode_slice(&values), [7, 4, 7, 5, 4]);
        assert_eq!(field.decode(lagrange_interpolation_at_zero(&points, &values, field)), 4);
    }

    macro_rules! all_tests {
        ($field:ty) => {
            use super::*;
            #[test] fn test_lagrange_interpolation_at_zero() { super::test_lagrange_interpolation_at_zero::<$field>(); }
        }
    }
    
    mod natural { all_tests!(::fields::NaturalPrimeField<i64>); }
    mod montgomery { all_tests!(::fields::MontgomeryField32); }
    #[cfg(feature="largefield")] mod large { all_tests!(::fields::LargePrimeField); }
    
}
