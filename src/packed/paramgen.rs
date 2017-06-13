
#![doc(hidden)]

//! Optional helper methods for parameter generation

extern crate num_traits;
extern crate primal;

use std::borrow::Borrow;
use std::ops::*;
use self::num_traits::{Zero, One};

#[cfg_attr(rustfmt, rustfmt_skip)]
fn check_prime_form<I, J: Borrow<I>>(min_p: J, n: J, m: J, p: J) -> bool 
where 
    I: PartialOrd + Zero + One,
    for<'l> &'l     I: Sub<I, Output=I>,
    for<'l, 'r> &'l I: Rem<&'r I, Output=I>,
                    I: Div<Output=I>,
    for<'l, 'r> &'l I: Mul<&'r I, Output=I>,
{
    if p.borrow() < min_p.borrow() { return false; }

    let q = p.borrow() - I::one();
    if &q % n.borrow() != I::zero() { return false; }
    if &q % m.borrow() != I::zero() { return false; }

    let k = q / (n.borrow() * m.borrow());
    if &k % n.borrow() == I::zero() { return false; }
    if &k % m.borrow() == I::zero() { return false; }

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
    assert_eq!(factor(41), vec![] as Vec<usize>);
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
