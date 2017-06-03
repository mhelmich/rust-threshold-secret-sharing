// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.

//! Packed (or ramp) variant of Shamir secret sharing,
//! allowing efficient sharing of several secrets together.

use rand;
use fields::{Field, Encode};

/// Parameters for the packed variant of Shamir secret sharing,
/// specifying number of secrets shared together, total number of shares, and privacy threshold.
///
/// This scheme generalises
/// [Shamir's scheme](https://en.wikipedia.org/wiki/Shamir%27s_Secret_Sharing)
/// by simultaneously sharing several secrets, at the expense of leaving a gap
/// between the privacy threshold and the reconstruction limit.
///
/// The Fast Fourier Transform is used for efficiency reasons,
/// allowing most operations run to quasilinear time `O(n.log(n))` in `share_count`.
/// An implication of this is that secrets and shares are positioned on positive powers of
/// respectively an `n`-th and `m`-th principal root of unity,
/// where `n` is a power of 2 and `m` a power of 3.
///
/// As a result there exist several constraints between the various parameters:
///
/// * `prime` must be a prime large enough to hold the secrets we plan to share
/// * `share_count` must be at least `secret_count + threshold` (the reconstruction limit)
/// * `secret_count + threshold + 1` must be a power of 2
/// * `share_count + 1` must be a power of 3
/// * `omega_secrets` must be a `(secret_count + threshold + 1)`-th root of unity
/// * `omega_shares` must be a `(share_count + 1)`-th root of unity
///
/// An optional `paramgen` feature provides methods for finding suitable parameters satisfying
/// these somewhat complex requirements, in addition to several fixed parameter choices.
#[derive(Debug,Clone,PartialEq)]
pub struct PackedSecretSharing<F: Field>
{
    // abstract properties

    /// Maximum number of shares that can be known without exposing the secrets
    /// (privacy threshold).
    pub threshold: usize,
    /// Number of shares to split the secrets into.
    pub share_count: usize,
    /// Number of secrets to share together.
    pub secret_count: usize,

    // implementation configuration

    /// Finite field in which computation is taking place.
    pub field: F,
    /// `m`-th principal root of unity in Zp, where `m = secret_count + threshold + 1`
    /// must be a power of 2.
    pub omega_secrets: F::E,
    /// `n`-th principal root of unity in Zp, where `n = share_count + 1` must be a power of 3.
    pub omega_shares: F::E,
}

impl<F> PackedSecretSharing<F> 
where F: Field, F: Encode<u32>, F::E: Clone
{
    /// Minimum number of shares required to reconstruct secrets.
    ///
    /// For this scheme this is always `secret_count + threshold`
    pub fn reconstruct_limit(&self) -> usize {
        self.threshold + self.secret_count
    }

    /// Generate `share_count` shares for the `secrets` vector.
    ///
    /// The length of `secrets` must be `secret_count`.
    /// It is safe to pad with anything, including zeros.
    pub fn share(&self, secrets: &[F::E]) -> Vec<F::E> {
        assert_eq!(secrets.len(), self.secret_count);
        // sample polynomial
        let mut poly = self.sample_polynomial(secrets);
        assert_eq!(poly.len(), self.reconstruct_limit() + 1);
        // .. and extend it (with zeroes)
        poly.extend(vec![self.field.zero(); self.share_count - self.reconstruct_limit()]);
        assert_eq!(poly.len(), self.share_count + 1);
        // evaluate polynomial to generate shares
        let mut shares = self.evaluate_polynomial(poly);
        // .. but remove first element since it should not be used as a share (it's always zero)
        assert!(self.field.eq(&shares[0], self.field.zero()));
        shares.remove(0);
        // return
        assert_eq!(shares.len(), self.share_count);
        shares
    }
    
    pub fn deterministic_share(&self, secrets_and_randomness: &[F::E]) -> Vec<F::E> {
        let mut values = secrets_and_randomness.to_vec();
        values.insert(0, self.field.zero());
        assert_eq!(values.len(), self.reconstruct_limit() + 1);
        ::numtheory::fft::fft2_inverse(&self.field, &mut *values, self.omega_secrets.clone());
        let mut poly = values;
        
        // TODO unify this with `share`
        
        assert_eq!(poly.len(), self.reconstruct_limit() + 1);
        // .. and extend it (with zeroes)
        poly.extend(vec![self.field.zero(); self.share_count - self.reconstruct_limit()]);
        assert_eq!(poly.len(), self.share_count + 1);
        // evaluate polynomial to generate shares
        let mut shares = self.evaluate_polynomial(poly);
        // .. but remove first element since it should not be used as a share (it's always zero)
        assert!(self.field.eq(&shares[0], self.field.zero()));
        shares.remove(0);
        // return
        assert_eq!(shares.len(), self.share_count);
        shares
    }

    fn sample_polynomial(&self, secrets: &[F::E]) -> Vec<F::E> {
        assert_eq!(secrets.len(), self.secret_count);
        // sample randomness using secure randomness
        let mut rng = rand::OsRng::new().unwrap();
        let randomness = self.field.sample_with_replacement(self.threshold, &mut rng);
        debug_assert!(self.field.neq(&randomness[0], &randomness[1])); // small probability for false negative
        // recover polynomial
        let coefficients = self.recover_polynomial(secrets, randomness);
        assert_eq!(coefficients.len(), self.reconstruct_limit() + 1);
        coefficients
    }

    fn recover_polynomial(&self, secrets: &[F::E], randomness: Vec<F::E>) -> Vec<F::E> {
        // fix the value corresponding to point 1 (zero)
        let mut values = vec![self.field.zero()];
        // let the subsequent values correspond to the secrets
        values.extend(secrets.to_vec()); // TODO can we do without to_vec?
        // fill in with random values
        values.extend(randomness);
        // run backward FFT to recover polynomial in coefficient representation
        assert_eq!(values.len(), self.reconstruct_limit() + 1);
        // in-place FFT to turn values into coefficients
        ::numtheory::fft::fft2_inverse(&self.field, &mut *values, self.omega_secrets.clone());
        let coefficients = values;
        coefficients
    }

    fn evaluate_polynomial(&self, mut coefficients: Vec<F::E>) -> Vec<F::E> {
        assert_eq!(coefficients.len(), self.share_count + 1);
        ::numtheory::fft::fft3(&self.field, &mut coefficients, self.omega_shares.clone());
        let points = coefficients;
        points
    }

    /// Reconstruct the secrets from a large enough subset of the shares.
    ///
    /// `indices` are the ranks of the known shares as output by the `share` method,
    ///  while `values` are the actual values of these shares.
    /// Both must have the same number of elements, and at least `reconstruct_limit`.
    ///
    /// The resulting vector is of length `secret_count`.
    pub fn reconstruct(&self, indices: &[u32], shares: &[F::E]) -> Vec<F::E> {
        assert!(shares.len() == indices.len());
        assert!(shares.len() >= self.reconstruct_limit());
        if shares.len() == self.share_count {
            // we're in the special case where we can use the FFTs for interpolation
            let mut values = shares.to_vec();
            values.insert(0, self.field.zero());
            ::numtheory::fft::fft3_inverse(&self.field, &mut values, self.omega_shares.clone());
            let mut coefficients = values.into_iter()
                .take(self.reconstruct_limit() + 1)
                .collect::<Vec<_>>();
            ::numtheory::fft::fft2(&self.field, &mut coefficients, self.omega_secrets.clone());
            let secrets = coefficients.into_iter()
                .skip(1)
                .take(self.secret_count)
                .collect();
            secrets
        } else {
            // we cannot use the FFT so default to Newton interpolation
            let mut points: Vec<F::E> = indices.iter()
                .map(|x| self.field.pow(&self.omega_shares, x + 1))
                .collect();
            let mut values = shares.to_vec();
            // insert missing value for point 1 (zero)
            points.insert(0, self.field.one());
            values.insert(0, self.field.zero());
            // interpolate using Newton's method
            // TODO optimise by using Newton-equally-space variant
            let poly = ::numtheory::newton_interpolation_general(&points, &values, &self.field);
            // evaluate at omega_secrets points to recover secrets
            // TODO optimise to avoid re-computation of power
            let secrets = (1..self.reconstruct_limit())
                .map(|e| self.field.pow(&self.omega_secrets, e as u32))
                .map(|point| ::numtheory::newton_evaluate(&poly, point, &self.field))
                .take(self.secret_count)
                .collect();
            secrets
        }
    }
    
    pub fn fully_reconstruct(&self, indices: &[u32], shares: &[F::E]) -> Vec<F::E> {
        // TODO unify code with `reconstruct` (only difference is how much is removed at end)
        
        assert!(shares.len() == indices.len());
        assert!(shares.len() >= self.reconstruct_limit());
        if shares.len() == self.share_count {
            // we're in the special case where we can use the FFTs for interpolation
            let mut values = shares.to_vec();
            values.insert(0, self.field.zero());
            ::numtheory::fft::fft3_inverse(&self.field, &mut values, self.omega_shares.clone());
            let mut coefficients = values.into_iter()
                .take(self.reconstruct_limit() + 1)
                .collect::<Vec<_>>();
            ::numtheory::fft::fft2(&self.field, &mut coefficients, self.omega_secrets.clone());
            let mut secrets = coefficients;
            secrets.remove(0);
            secrets
        } else {
            // we cannot use the FFT so default to Newton interpolation
            let mut points: Vec<F::E> = indices.iter()
                .map(|x| self.field.pow(&self.omega_shares, x + 1))
                .collect();
            let mut values = shares.to_vec();
            // insert missing value for point 1 (zero)
            points.insert(0, self.field.one());
            values.insert(0, self.field.zero());
            // interpolate using Newton's method
            // TODO optimise by using Newton-equally-space variant
            let poly = ::numtheory::newton_interpolation_general(&points, &values, &self.field);
            // evaluate at omega_secrets points to recover secrets
            // TODO optimise to avoid re-computation of power
            let secrets = (1..self.reconstruct_limit())
                .map(|e| self.field.pow(&self.omega_secrets, e as u32))
                .map(|point| ::numtheory::newton_evaluate(&poly, point, &self.field))
                .collect();
            secrets
        }
    }
}


mod instances 
{    
    use super::*;
    use fields::NaturalPrimeField;

    /// Example of tiny PSS settings, for sharing 3 secrets into 8 shares, with
    /// a privacy threshold of 4.
    pub static PSS_4_8_3: PackedSecretSharing<NaturalPrimeField<i64>> = PackedSecretSharing {
        threshold: 4,
        share_count: 8,
        secret_count: 3,
        field: NaturalPrimeField(433), // TODO
        omega_secrets: 354,
        omega_shares: 150,
    };

    /// Example of small PSS settings, for sharing 3 secrets into 26 shares, with
    /// a privacy threshold of 4.
    pub static PSS_4_26_3: PackedSecretSharing<NaturalPrimeField<i64>> = PackedSecretSharing {
        threshold: 4,
        share_count: 26,
        secret_count: 3,
        field: NaturalPrimeField(433), // TODO
        omega_secrets: 354,
        omega_shares: 17,
    };

    /// Example of PSS settings, for sharing 100 secrets into 728 shares, with
    /// a privacy threshold of 155.
    pub static PSS_155_728_100: PackedSecretSharing<NaturalPrimeField<i64>> = PackedSecretSharing {
        threshold: 155,
        share_count: 728,
        secret_count: 100,
        field: NaturalPrimeField(746497), // TODO
        omega_secrets: 95660,
        omega_shares: 610121,
    };

    /// Example of PSS settings, for sharing 100 secrets into 19682 shares, with
    /// a privacy threshold of 155.
    pub static PSS_155_19682_100: PackedSecretSharing<NaturalPrimeField<i64>> = PackedSecretSharing {
        threshold: 155,
        share_count: 19682,
        secret_count: 100,
        field: NaturalPrimeField(5038849), // TODO
        omega_secrets: 4318906,
        omega_shares: 1814687,
    };
}
pub use self::instances::*;


#[cfg(test)]
mod tests {

    use super::*;
    use ::fields::*;
    
    pub fn test_recover_polynomial<F>()
    where F: PrimeField + Encode<u32> + Decode<u32> + Clone, F::P: From<u32>, F::E: Clone
    {
        let field = F::new(433.into());
        let pss = PackedSecretSharing {
            threshold: 4,
            share_count: 8,
            secret_count: 3,
            omega_secrets: field.encode(354),
            omega_shares: field.encode(150),
            field: field.clone(),
        };
        
        let secrets = vec![1, 2, 3];
        let randomness = vec![8, 8, 8, 8];  // use fixed randomness
        let poly = pss.recover_polynomial(&field.encode_slice(secrets), field.encode_slice(randomness));
        assert_eq!(
            field.decode_slice(poly),
            vec![113, 51, 261, 267, 108, 432, 388, 112]
        );
    }
    
    #[cfg_attr(rustfmt, rustfmt_skip)]
    pub fn test_evaluate_polynomial<F>()
    where F: PrimeField + Encode<u32> + Decode<u32> + Clone, F::P: From<u32>, F::E: Clone
    {
        let field = F::new(433.into());
        let pss = PackedSecretSharing {
            threshold: 4,
            share_count: 26,
            secret_count: 3,
            omega_secrets: field.encode(354),
            omega_shares: field.encode(17),
            field: field.clone(),
        };
        
        let poly = field.encode_slice([113,  51, 261, 267, 108, 432, 388, 112,   0,
                                         0,   0,   0,   0,   0,   0,   0,   0,   0,
                                         0,   0,   0,   0,   0,   0,   0,   0,   0]);
        let points = &pss.evaluate_polynomial(poly);
        assert_eq!(
            field.decode_slice(points),
            [   0, 77, 230,  91, 286, 179, 337,  83, 212,
               88, 406, 58, 425, 345, 350, 336, 430, 404,
               51, 60, 305, 395,  84, 156, 160, 112, 422]
        );
    }
    
}

macro_rules! all_packed_tests {
    ($field:ty) => {
        #[test] fn test_recover_polynomial() { ::packed::tests::test_recover_polynomial::<$field>(); }
        #[test] fn test_evaluate_polynomial() { ::packed::tests::test_evaluate_polynomial::<$field>(); }
    }
}

#[cfg(test)] mod natural    { all_packed_tests!(::fields::NaturalPrimeField<i64>); }
#[cfg(test)] mod montgomery { all_packed_tests!(::fields::MontgomeryField32); }
#[cfg(all(test, feature="largefield"))] mod large { all_packed_tests!(::fields::LargePrimeField); }


#[cfg(test)]
mod old_tests {

    use super::*;
    use ::fields::*;

    #[test]
    #[cfg_attr(rustfmt, rustfmt_skip)]
    fn test_share() {
        let ref pss = PSS_4_26_3;
        let ref field = pss.field;

        // do sharing
        let secrets = field.encode_slice([5, 6, 7]);
        let mut shares = pss.share(&secrets);

        // manually recover secrets
        use numtheory::mod_evaluate_polynomial;
        shares.insert(0, 0);
        ::numtheory::fft::fft3_inverse(field, &mut *shares, pss.omega_shares);
        let poly = shares;
        let recovered_secrets: Vec<i64> = (1..secrets.len() + 1)
            .map(|i| pss.field.pow(field.encode(pss.omega_secrets as u32), i as u32))
            .map(|point| mod_evaluate_polynomial(&poly, point, field))
            .collect();

        assert_eq!(
            field.decode_slice(recovered_secrets),
            field.decode_slice(secrets)
        );
    }

    #[test]
    fn test_large_share() {
        let ref pss = PSS_155_19682_100;
        let secrets = vec![5 ; pss.secret_count];
        let shares = pss.share(&secrets);
        assert_eq!(shares.len(), pss.share_count);
    }

    #[test]
    fn test_share_reconstruct() {
        let ref pss = PSS_4_26_3;
        let secrets = vec![5, 6, 7];
        let shares = pss.share(&pss.field.encode_slice(&secrets));

        // reconstruction must work for all shares
        let indices: Vec<u32> = (0..shares.len() as u32).collect();
        let recovered_secrets = pss.reconstruct(&indices, &shares);
        assert_eq!(pss.field.decode_slice(recovered_secrets), secrets);

        // .. and for only sufficient shares
        let indices: Vec<u32> = (0..pss.reconstruct_limit() as u32).collect();
        let recovered_secrets = pss.reconstruct(&indices, &shares[0..pss.reconstruct_limit()]);
        print!("lenght is {:?}", indices.len());
        assert_eq!(pss.field.decode_slice(recovered_secrets), secrets);
    }

    #[test]
    fn test_share_additive_homomorphism() {
        let ref pss = PSS_4_26_3;

        let secrets_1 = vec![1, 2, 3];
        let secrets_2 = vec![4, 5, 6];
        let shares_1 = pss.share(&secrets_1);
        let shares_2 = pss.share(&secrets_2);

        // add shares pointwise
        let shares_sum: Vec<i64> =
            shares_1.iter().zip(shares_2).map(|(a, b)| (a + b) % pss.field.0).collect(); // TODO

        // reconstruct sum, using same reconstruction limit
        let reconstruct_limit = pss.reconstruct_limit();
        let indices: Vec<u32> = (0..reconstruct_limit as u32).collect();
        let shares = &shares_sum[0..reconstruct_limit];
        let recovered_secrets = pss.reconstruct(&indices, shares);

        assert_eq!(pss.field.decode_slice(recovered_secrets), [5, 7, 9]);
    }

    #[test]
    fn test_share_multiplicative_homomorphism() {
        let ref pss = PSS_4_26_3;

        let secrets_1 = vec![1, 2, 3];
        let secrets_2 = vec![4, 5, 6];
        let shares_1 = pss.share(&secrets_1);
        let shares_2 = pss.share(&secrets_2);

        // multiply shares pointwise
        let shares_product: Vec<i64> =
            shares_1.iter().zip(shares_2).map(|(a, b)| (a * b) % pss.field.0).collect(); // TODO

        // reconstruct product, using double reconstruction limit
        let reconstruct_limit = pss.reconstruct_limit() * 2;
        let indices: Vec<u32> = (0..reconstruct_limit as u32).collect();
        let shares = &shares_product[0..reconstruct_limit];
        let recovered_secrets = pss.reconstruct(&indices, shares);

        assert_eq!(pss.field.decode_slice(recovered_secrets), [4, 10, 18]);
    }

}


#[cfg(feature = "paramgen")]
mod paramgen;
#[cfg(feature = "paramgen")]
pub use self::paramgen::*;
