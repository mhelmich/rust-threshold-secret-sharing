// Copyright (c) 2017 rust-threshold-secret-sharing developers

#![allow(dead_code)]

#[macro_use]
extern crate bencher;
extern crate threshold_secret_sharing as tss;

use tss::*;
use bencher::Bencher;

pub trait Config {
    fn privacy() -> usize;
    fn shares() -> usize;
    fn prime() -> u32;
    fn omega() -> u32;
}

struct Tiny;
impl Config for Tiny {
    fn privacy() -> usize { 1 }
    fn shares() -> usize { 2 }
    fn prime() -> u32 { 433 }
    fn omega() -> u32 { 198 }
}

struct Small;
impl Config for Small {
    fn privacy() -> usize { 2 }
    fn shares() -> usize { 8 }
    fn prime() -> u32 { 433 }
    fn omega() -> u32 { 150 }
}

struct Medium;
impl Config for Medium {
    fn privacy() -> usize { 7 }
    fn shares() -> usize { 26 }
    fn prime() -> u32 { 433 }
    fn omega() -> u32 { 17 }
}

struct Large;
impl Config for Large {
    fn privacy() -> usize { 20 }
    fn shares() -> usize { 80 }
    fn prime() -> u32 { 746497 }
    fn omega() -> u32 { 69177 }
}

struct Huge;
impl Config for Huge {
    fn privacy() -> usize { 61 }
    fn shares() -> usize { 242 }
    fn prime() -> u32 { 746497 }
    fn omega() -> u32 { 595577 }
}

pub fn share_fft<C: Config, F>(b: &mut Bencher) 
where F: PrimeField + New<u32> + Encode<u32>, F::P: From<u32>, F::E: Clone
{
    let ref field = F::new(C::prime());
    let ref omega_shares = field.encode(C::omega());
    let ref values = field.encode_slice(vec![5; C::privacy() + 1]);
    
    b.iter(|| {
        let mut data = values.clone();
        data.extend(vec![field.zero(); C::shares() - C::privacy()]);
        ::numtheory::fft::fft3(field, &mut *data, omega_shares);
        let _shares = data;
    });
}

pub fn share_horner<C: Config, F>(b: &mut Bencher) 
where F: PrimeField + New<u32> + Encode<u32>, F::P: From<u32>, F::E: Clone
{
    let ref field = F::new(C::prime());
    let ref values = field.encode_slice(vec![5; C::privacy() + 1]);
    
    b.iter(|| {
        let _shares = (1..C::shares() as u32+1)
            .map(|p| ::numtheory::mod_evaluate_polynomial(values, field.encode(p), field))
            .collect::<Vec<_>>();
    });
}

benchmark_group!(group
    , share_fft    <Tiny, MontgomeryField32>
    , share_horner <Tiny, MontgomeryField32>
    , share_fft    <Small, MontgomeryField32>
    , share_horner <Small, MontgomeryField32>
    , share_fft    <Medium, MontgomeryField32>
    , share_horner <Medium, MontgomeryField32>
    , share_fft    <Large, MontgomeryField32>
    , share_horner <Large, MontgomeryField32>
    , share_fft    <Huge, MontgomeryField32>
    , share_horner <Huge, MontgomeryField32>
);

benchmark_main!(group);
