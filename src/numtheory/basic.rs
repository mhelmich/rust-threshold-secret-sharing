// Copyright (c) 2017 rust-threshold-secret-sharing developers

use fields::Field;
use std::borrow::Borrow;

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

pub fn generic_mod_pow<F>(field: &F, a: F::E, e: u32) -> F::E
where
    F: Field,
{
    // TODO improve (or at least compare to non-generic)

    let mut x = a;
    let mut e = e;
    let mut acc = field.one();
    while e > 0 {
        if e % 2 == 0 {
            // even
            // no-op
        } else {
            // odd
            acc = field.mul(&acc, &x);
        }
        x = field.mul(&x, &x); // waste one of these by having it here but code is simpler (tiny bit)
        e = e >> 1;
    }
    acc
}

/// Evaluate polynomial given by `coefficients` at `point`.
pub fn mod_evaluate_polynomial<F, E>(coefficients: &[F::E], point: E, field: &F) -> F::E
where
    F: Field,
    E: Borrow<F::E>,
{
    // evaluate using Horner's rule
    //  - to combine with fold we consider the coefficients in reverse order
    coefficients
        .iter()
        .rev()
        .fold(field.zero(), |partial, coef| {
            field.add(field.mul(partial, point.borrow()), coef)
        })
}

pub fn weighted_sum<F>(values: &[F::E], weights: &[F::E], field: &F) -> F::E
where
    F: Field,
{
    values
        .iter()
        .zip(weights)
        .map(|(v, w)| field.mul(v, w))
        .fold(field.zero(), |sum, term| field.add(sum, term))
}

#[cfg(test)]
mod tests {

    use super::*;
    use fields;

    #[test]
    fn test_mod_pow() {
        assert_eq!(mod_pow(2, 0, 17), 1);
        assert_eq!(mod_pow(2, 3, 17), 8);
        assert_eq!(mod_pow(2, 6, 17), 13);

        assert_eq!(mod_pow(-3, 0, 17), 1);
        assert_eq!(mod_pow(-3, 1, 17), -3);
        assert_eq!(mod_pow(-3, 15, 17), -6);
    }

    #[test]
    fn test_mod_evaluate_polynomial() {
        let poly = vec![1, 2, 3, 4, 5, 6];
        let point = 5;
        let field = fields::NaturalPrimeField(17);
        assert_eq!(mod_evaluate_polynomial(&poly, point, &field), 4);
    }
}
