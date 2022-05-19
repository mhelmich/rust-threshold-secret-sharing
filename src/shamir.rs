// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.

//! Standard [Shamir secret sharing](https://en.wikipedia.org/wiki/Shamir%27s_Secret_Sharing)
//! for a single secret.

use rand;

use fields::Encode;
use fields::Field;

/// Parameters for the Shamir scheme, specifying privacy threshold and total number of shares.
///
/// There are very few constraints except for the obvious ones:
///
/// * `prime` must be a prime large enough to hold the secrets we plan to share
/// * `share_count` must be at least `threshold + 1` (the reconstruction limit)
///
/// # Example:
///
/// ```
///    use threshold_secret_sharing::*;
///    let tss = ShamirSecretSharing {
///        threshold: 9,
///        share_count: 20,
///        field: NaturalPrimeField(41)
///    };
///
///    let secret = 5;
///    let all_shares = tss.share(secret);
///
///    let reconstruct_share_count = tss.reconstruct_limit();
///
///    let indices: Vec<usize> = (0..reconstruct_share_count).collect();
///    let shares: &[i64] = &all_shares[0..reconstruct_share_count];
///    let recovered_secret = tss.reconstruct(&indices, shares);
///
///    println!("The recovered secret is {}", recovered_secret);
///    assert_eq!(recovered_secret, secret);
/// ```
#[derive(Debug)]
pub struct ShamirSecretSharing<F>
where
    F: Field,
    F::E: Clone,
{
    /// Maximum number of shares that can be known without exposing the secret.
    pub threshold: usize,
    /// Number of shares to split the secret into.
    pub share_count: usize,
    /// Finite field in which computation takes place.
    pub field: F,
}

impl<F> ShamirSecretSharing<F>
where
    F: Field,
    F: Encode<u32>,
    F::E: Clone,
{
    /// Minimum number of shares required to reconstruct secret.
    ///
    /// For this scheme this is always `threshold + 1`.
    pub fn reconstruct_limit(&self) -> usize {
        self.threshold + 1
    }

    /// Generate `share_count` shares from `secret`.
    pub fn share(&self, secret: F::E) -> Vec<F::E> {
        let poly = self.sample_polynomial(secret);
        self.evaluate_polynomial(&poly)
    }

    fn sample_polynomial(&self, zero_value: F::E) -> Vec<F::E> {
        // fix the first coefficient (corresponding to the evaluation at zero)
        let mut coefficients = vec![zero_value];
        // sample the remaining coefficients randomly using secure randomness
        let mut rng = rand::OsRng::new().unwrap();
        let random_coefficients = self.field.sample_with_replacement(self.threshold, &mut rng);
        coefficients.extend(random_coefficients);
        // return
        coefficients
    }

    fn evaluate_polynomial(&self, coefficients: &[F::E]) -> Vec<F::E> {
        // evaluate at all points
        (1..self.share_count + 1)
            .map(|point| {
                ::numtheory::mod_evaluate_polynomial(
                    coefficients,
                    self.field.encode(point as u32),
                    &self.field,
                )
            })
            .collect()
    }

    /// Reconstruct `secret` from a large enough subset of the shares.
    ///
    /// `indices` are the ranks of the known shares as output by the `share` method,
    /// while `values` are the actual values of these shares.
    /// Both must have the same number of elements, and at least `reconstruct_limit`.
    pub fn reconstruct(&self, indices: &[usize], shares: &[F::E]) -> F::E {
        assert!(shares.len() == indices.len());
        assert!(shares.len() >= self.reconstruct_limit());
        // add one to indices to get points
        let points: Vec<F::E> = indices
            .iter()
            .map(|&i| {
                self.field
                    .add(self.field.encode(i as u32), self.field.one())
            })
            .collect();
        // interpolate
        ::numtheory::lagrange_interpolation_at_zero(&*points, &shares, &self.field)
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use fields::NaturalPrimeField;

    // Small preset parameters for tests.
    pub static SHAMIR_5_20: ShamirSecretSharing<NaturalPrimeField<i64>> = ShamirSecretSharing {
        threshold: 5,
        share_count: 20,
        field: NaturalPrimeField(41),
    };

    #[test]
    fn test_evaluate_polynomial() {
        let ref tss = SHAMIR_5_20;
        let poly = vec![1, 2, 0];
        let values = tss.evaluate_polynomial(&poly);
        assert_eq!(
            *values,
            [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 0]
        );
    }

    #[test]
    fn wikipedia_example() {
        let tss = ShamirSecretSharing {
            threshold: 2,
            share_count: 6,
            field: NaturalPrimeField(1613),
        };

        let shares = tss.evaluate_polynomial(&[1234, 166, 94]);
        assert_eq!(&*shares, &[1494, 329, 965, 176, 1188, 775]);

        assert_eq!(tss.reconstruct(&[0, 1, 2], &shares[0..3]), 1234);
        assert_eq!(tss.reconstruct(&[1, 2, 3], &shares[1..4]), 1234);
        assert_eq!(tss.reconstruct(&[2, 3, 4], &shares[2..5]), 1234);
    }

    #[test]
    fn test_shamir() {
        let tss = ShamirSecretSharing {
            threshold: 2,
            share_count: 6,
            field: NaturalPrimeField(41),
        };
        let secret = 1;
        let shares = tss.share(secret);
        assert_eq!(tss.reconstruct(&[0, 1, 2], &shares[0..3]), secret);
        assert_eq!(tss.reconstruct(&[1, 2, 3], &shares[1..4]), secret);
        assert_eq!(tss.reconstruct(&[2, 3, 4, 5], &shares[2..6]), secret);
    }
}
