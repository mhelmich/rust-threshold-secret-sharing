// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.

//! Various number theoretic utility functions used in the library.

pub mod fft;

pub mod newton;
pub use self::newton::*;

pub mod lagrange;
pub use self::lagrange::*;

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
where F: Field
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
        x = field.mul(&x, &x);  // waste one of these by having it here but code is simpler (tiny bit)
        e = e >> 1;
    }
    acc
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
where F: Field, F: Encode<u32>, F::E: Clone
{
    // evaluate using Horner's rule
    //  - to combine with fold we consider the coefficients in reverse order
    coefficients.iter().rev()
        .fold(field.zero(), |partial, coef| field.add(field.mul(partial, &point), coef))
}


#[cfg(test)]
mod tests {
    
    use super::*;
    use fields;
    
    #[test]
    fn test_gcd() {
        assert_eq!(gcd(12, 16), (4, -1, 1));
    }
    
    #[test]
    pub fn test_binary_egcd() {
        assert_eq!(binary_egcd(10, 4), (2, 1, -2));
    }
    
    #[test]
    fn test_mod_inverse() {
        assert_eq!(mod_inverse(3, 7), 5);
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
    
    #[test]
    fn test_mod_evaluate_polynomial() {
        let poly = vec![1, 2, 3, 4, 5, 6];
        let point = 5;
        let field = fields::NaturalPrimeField(17);
        assert_eq!(mod_evaluate_polynomial(&poly, point, &field), 4);
    }
    
}
