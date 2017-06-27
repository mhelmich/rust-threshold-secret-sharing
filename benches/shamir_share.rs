// Copyright (c) 2017 rust-threshold-secret-sharing developers

#![allow(dead_code)]

#[macro_use]
extern crate bencher;
extern crate threshold_secret_sharing as tss;

use tss::*;
use bencher::Bencher;

pub trait Config {
    fn prime() -> u32;
    fn count_secret_randomness() -> usize;
    fn count_shares() -> usize;
    fn omega_shares() -> u32;
}

struct Tiny;
impl Config for Tiny {
    fn prime() -> u32 { 433 }
    fn count_secret_randomness() -> usize { 2 + 1 }
    fn count_shares() -> usize { 3 }
    fn omega_shares() -> u32 { 198 }
}

struct Small;
impl Config for Small {
    fn prime() -> u32 { 433 }
    fn count_secret_randomness() -> usize { 5 + 1 }
    fn count_shares() -> usize { 9 }
    fn omega_shares() -> u32 { 150 }
}

struct Medium;
impl Config for Medium {
    fn prime() -> u32 { 746497 }
    fn count_secret_randomness() -> usize { 365 + 1 }
    fn count_shares() -> usize { 729 }
    fn omega_shares() -> u32 { 610121 }
}

struct Large;
impl Config for Large {
    fn prime() -> u32 { 416544769 }
    fn count_secret_randomness() -> usize { 1094 + 1 }
    fn count_shares() -> usize { 2187 }
    fn omega_shares() -> u32 { 35260261 }
}

pub fn share_fft<C: Config, F>(b: &mut Bencher) 
where F: PrimeField + New<u32> + Encode<u32>, F::P: From<u32>, F::E: Clone
{
    let ref field = F::new(C::prime());
    let ref omega_shares = field.encode(C::omega_shares());
    let ref values = field.encode_slice(vec![5; C::count_secret_randomness()]);
    
    b.iter(|| {
        let mut data = values.clone();
        data.extend(vec![field.zero(); C::count_shares() - C::count_secret_randomness()]);
        ::numtheory::fft::fft3(field, &mut *data, omega_shares);
        let _shares = data;
    });
}

pub fn share_horner<C: Config, F>(b: &mut Bencher) 
where F: PrimeField + New<u32> + Encode<u32>, F::P: From<u32>, F::E: Clone
{
    let ref field = F::new(C::prime());
    let ref values = field.encode_slice(vec![5; C::count_secret_randomness()]);
    
    b.iter(|| {
        let _shares = (1..C::count_shares() as u32+1)
            .map(|p| ::numtheory::mod_evaluate_polynomial(values, field.encode(p), field))
            .collect::<Vec<_>>();
    });
}

benchmark_group!(tiny
    , share_fft    <Tiny, MontgomeryField32>
    , share_horner <Tiny, MontgomeryField32>
);

benchmark_group!(small
    , share_fft    <Small, MontgomeryField32>
    , share_horner <Small, MontgomeryField32>
);
 
benchmark_group!(medium
    , share_fft    <Medium, MontgomeryField32>
    , share_horner <Medium, MontgomeryField32>
);

benchmark_group!(large
    , share_fft    <Large, MontgomeryField32>
    , share_horner <Large, MontgomeryField32>
);

benchmark_main!(
    tiny
    , small
    , medium
    , large
);
