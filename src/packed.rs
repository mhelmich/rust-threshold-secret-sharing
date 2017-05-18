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
#[derive(Debug,Copy,Clone,PartialEq)]
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
where F: Field, F: Encode<u32>, F::E: Copy
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
        assert!(self.field.eq(shares[0], self.field.zero()));
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
        debug_assert!(self.field.neq(randomness[0], randomness[1])); // small probability for false negative
        // recover polynomial
        let coefficients = self.recover_polynomial(secrets, randomness);
        assert_eq!(coefficients.len(), self.reconstruct_limit() + 1);
        coefficients
    }

    fn recover_polynomial(&self, secrets: &[F::E], randomness: Vec<F::E>) -> Vec<F::E> {
        // fix the value corresponding to point 1 (zero)
        let mut values = vec![self.field.zero()];
        // let the subsequent values correspond to the secrets
        values.extend(secrets);
        // fill in with random values
        values.extend(randomness);
        // run backward FFT to recover polynomial in coefficient representation
        assert_eq!(values.len(), self.reconstruct_limit() + 1);
        // in-place FFT to turn values into coefficients
        ::numtheory::fft::fft2_inverse(&self.field, &mut *values, self.omega_secrets);
        let coefficients = values;
        coefficients
    }

    fn evaluate_polynomial(&self, mut coefficients: Vec<F::E>) -> Vec<F::E> {
        assert_eq!(coefficients.len(), self.share_count + 1);
        ::numtheory::fft::fft3(&self.field, &mut coefficients, self.omega_shares);
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
        let mut points: Vec<F::E> = indices.iter()
            .map(|&x| self.field.pow(self.omega_shares, x + 1))
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
            .map(|e| self.field.pow(self.omega_secrets, e as u32))
            .map(|point| ::numtheory::newton_evaluate(&poly, point, &self.field))
            .take(self.secret_count)
            .collect();
        secrets
    }
}


mod instances {
    
    use super::*;
    use fields::natural::NaturalPrimeField;

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


#[doc(hidden)]
#[cfg(feature = "paramgen")]
pub mod paramgen {

    //! Optional helper methods for parameter generation

    extern crate primal;

    #[cfg_attr(rustfmt, rustfmt_skip)]
    fn check_prime_form(min_p: usize, n: usize, m: usize, p: usize) -> bool {
        if p < min_p { return false; }

        let q = p - 1;
        if q % n != 0 { return false; }
        if q % m != 0 { return false; }

        let q = q / (n * m);
        if q % n == 0 { return false; }
        if q % m == 0 { return false; }

        return true;
    }

    #[test]
    fn test_check_prime_form() {
        assert_eq!(primal::Primes::all().find(|p| check_prime_form(198, 8, 9, *p)).unwrap(), 433);
    }

    fn factor(p: usize) -> Vec<usize> {
        let mut factors = vec![];
        let bound = (p as f64).sqrt().ceil() as usize;
        for f in 2..bound + 1 {
            if p % f == 0 {
                factors.push(f);
                factors.push(p / f);
            }
        }
        factors
    }

    #[test]
    fn test_factor() {
        assert_eq!(factor(40), [2, 20, 4, 10, 5, 8]);
        assert_eq!(factor(41), []);
    }

    fn find_field(min_p: usize, n: usize, m: usize) -> Option<(i64, i64)> {
        // find prime of right form
        let p = primal::Primes::all().find(|p| check_prime_form(min_p, n, m, *p)).unwrap();
        // find (any) generator
        let factors = factor(p - 1);
        for g in 2..p {
            // test generator against all factors of p-1
            let is_generator = factors.iter().all(|f| {
                let e = (p - 1) / f;
                ::numtheory::mod_pow(g as i64, e as u32, p as i64) != 1  // TODO check for negative value
            });
            // return
            if is_generator {
                return Some((p as i64, g as i64));
            }
        }
        // didn't find any
        None
    }

    #[test]
    fn test_find_field() {
        assert_eq!(find_field(198, 2usize.pow(3), 3usize.pow(2)).unwrap(), (433, 5));
        assert_eq!(find_field(198, 2usize.pow(3), 3usize.pow(3)).unwrap(), (433, 5));
        assert_eq!(find_field(198, 2usize.pow(8), 3usize.pow(6)).unwrap(), (746497, 5));
        assert_eq!(find_field(198, 2usize.pow(8), 3usize.pow(9)).unwrap(), (5038849, 29));

        // assert_eq!(find_field(198, 2usize.pow(11), 3usize.pow(8)).unwrap(), (120932353, 5));
        // assert_eq!(find_field(198, 2usize.pow(13), 3usize.pow(9)).unwrap(), (483729409, 23));
    }

    fn find_roots(n: usize, m: usize, p: i64, g: i64) -> (i64, i64) {
        let omega_secrets = ::numtheory::mod_pow(g, ((p - 1) / n as i64) as u32, p);
        let omega_shares = ::numtheory::mod_pow(g, ((p - 1) / m as i64) as u32, p);
        (omega_secrets, omega_shares)
    }

    #[test]
    fn test_find_roots() {
        assert_eq!(find_roots(2usize.pow(3), 3usize.pow(2), 433, 5), (354, 150));
        assert_eq!(find_roots(2usize.pow(3), 3usize.pow(3), 433, 5), (354, 17));
    }

    #[doc(hidden)]
    pub fn generate_parameters(min_size: usize, n: usize, m: usize) -> (i64, i64, i64) {
        // TODO settle option business once and for all (don't remember it as needed)
        let (prime, g) = find_field(min_size, n, m).unwrap();
        let (omega_secrets, omega_shares) = find_roots(n, m, prime, g);
        (prime, omega_secrets, omega_shares)
    }

    #[test]
    fn test_generate_parameters() {
        assert_eq!(generate_parameters(200, 2usize.pow(3), 3usize.pow(2)), (433, 354, 150));
        assert_eq!(generate_parameters(200, 2usize.pow(3), 3usize.pow(3)), (433, 354, 17));
    }

    fn is_power_of(x: usize, e: usize) -> bool {
        let power = (x as f64).log(e as f64).floor() as u32;
        e.pow(power) == x
    }

    #[test]
    fn test_is_power_of() {
        assert_eq!(is_power_of(4, 2), true);
        assert_eq!(is_power_of(5, 2), false);
        assert_eq!(is_power_of(6, 2), false);
        assert_eq!(is_power_of(7, 2), false);
        assert_eq!(is_power_of(8, 2), true);

        assert_eq!(is_power_of(4, 3), false);
        assert_eq!(is_power_of(5, 3), false);
        assert_eq!(is_power_of(6, 3), false);
        assert_eq!(is_power_of(7, 3), false);
        assert_eq!(is_power_of(8, 3), false);
        assert_eq!(is_power_of(9, 3), true);
    }
    
    use super::*;
    use ::fields::PrimeField;

    impl<F> PackedSecretSharing<F> 
    where F: PrimeField, F: Encode<u32>, F::P: From<u32>
    {

        /// Find suitable parameters with as small a prime field as possible.
        pub fn new(threshold: usize,
                   secret_count: usize,
                   share_count: usize)
                   -> PackedSecretSharing<F> {
            let min_size = share_count + secret_count + threshold + 1;
            Self::new_with_min_size(threshold, secret_count, share_count, min_size)
        }

        /// Find suitable parameters with a prime field of at least the specified size.
        pub fn new_with_min_size(threshold: usize,
                                 secret_count: usize,
                                 share_count: usize,
                                 min_size: usize)
                                 -> PackedSecretSharing<F> {

            let m = threshold + secret_count + 1;
            let n = share_count + 1;
            assert!(is_power_of(m, 2));
            assert!(is_power_of(n, 3));
            assert!(min_size >= share_count + secret_count + threshold + 1);

            let (prime, omega_secrets, omega_shares) = generate_parameters(min_size, m, n);
            
            let field = F::new((prime as u32).into());
            PackedSecretSharing {
                threshold: threshold,
                share_count: share_count,
                secret_count: secret_count,
                omega_secrets: field.encode(omega_secrets as u32),
                omega_shares: field.encode(omega_shares as u32),
                field: field,
            }
        }
    }

    #[test]
    fn test_new() {
        assert_eq!(PackedSecretSharing::new(155, 100, 728), super::PSS_155_728_100);
        assert_eq!(PackedSecretSharing::new_with_min_size(4, 3, 8, 200), super::PSS_4_8_3);
        assert_eq!(PackedSecretSharing::new_with_min_size(4, 3, 26, 200), super::PSS_4_26_3);
    }

}


#[cfg(test)]
mod tests {

    use super::*;
    use numtheory::*;
    use ::fields::Encode;
    
    #[test]
    fn test_recover_polynomial() {
        let ref pss = PSS_4_8_3;
        let secrets = vec![1, 2, 3];
        let randomness = vec![8, 8, 8, 8];  // use fixed randomness
        let poly = pss.recover_polynomial(&secrets, randomness);
        assert_eq!(
            positivise(&poly, pss.field.0),
            positivise(&[113, -382, -172, 267, -325, 432, 388, -321], pss.field.0)
        );
    }

    #[test]
    #[cfg_attr(rustfmt, rustfmt_skip)]
    fn test_evaluate_polynomial() {
        let ref pss = PSS_4_26_3;
        let poly = vec![113,  51, 261, 267, 108, 432, 388, 112,   0,
                          0,   0,   0,   0,   0,   0,   0,   0,   0,
                          0,   0,   0,   0,   0,   0,   0,   0,   0];
        let points = &pss.evaluate_polynomial(poly);
        assert_eq!(
            positivise(points, pss.field.0),
            vec![   0, 77, 230,  91, 286, 179, 337,  83, 212,
                   88, 406, 58, 425, 345, 350, 336, 430, 404,
                   51, 60, 305, 395,  84, 156, 160, 112, 422]
        );
    }

    #[test]
    #[cfg_attr(rustfmt, rustfmt_skip)]
    fn test_share() {
        let ref pss = PSS_4_26_3;
        let ref field = pss.field;

        // do sharing
        let secrets: Vec<i64> = vec![5, 6, 7].iter()
            .map(|&v| field.encode(v as u32))
            .collect();
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

        use numtheory::positivise;
        // TODO ideally eq for wrapping type should handle this, so no need to positivise
        assert_eq!(
            positivise(&recovered_secrets, field.0),
            secrets
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
        let shares = pss.share(&secrets);

        // reconstruction must work for all shares
        let indices: Vec<u32> = (0..shares.len() as u32).collect();
        let recovered_secrets = pss.reconstruct(&indices, &shares);
        assert_eq!(positivise(&recovered_secrets, pss.field.0), secrets);

        // .. and for only sufficient shares
        let indices: Vec<u32> = (0..pss.reconstruct_limit() as u32).collect();
        let recovered_secrets = pss.reconstruct(&indices, &shares[0..pss.reconstruct_limit()]);
        print!("lenght is {:?}", indices.len());
        assert_eq!(positivise(&recovered_secrets, pss.field.0), secrets);
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

        use numtheory::positivise;
        assert_eq!(positivise(&recovered_secrets, pss.field.0), vec![5, 7, 9]);
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

        use numtheory::positivise;
        assert_eq!(positivise(&recovered_secrets, pss.field.0), vec![4, 10, 18]);
    }

}
