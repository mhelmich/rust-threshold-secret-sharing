// Copyright (c) 2017 rust-threshold-secret-sharing developers

#![allow(dead_code)]

#[macro_use]
extern crate bencher;
extern crate threshold_secret_sharing as tss;

use tss::*;
use bencher::Bencher;

pub trait Config {
    fn prime() -> u32;
    fn count_secrets() -> usize;
    fn count_shares() -> usize;
    fn omega_secrets() -> u32;
    fn omega_shares() -> u32;
}

struct Small;
impl Config for Small {
    fn prime() -> u32 { 433 }
    fn count_secrets() -> usize { 8 }
    fn count_shares() -> usize { 9 }
    fn omega_secrets() -> u32 { 354 }
    fn omega_shares() -> u32 { 150 }
}

struct Medium;
impl Config for Medium {
    fn prime() -> u32 { 746497 }
    fn count_secrets() -> usize { 256 }
    fn count_shares() -> usize { 729 }
    fn omega_secrets() -> u32 { 95660 }
    fn omega_shares() -> u32 { 610121 }
}

struct Large;
impl Config for Large {
    fn prime() -> u32 { 416544769 }
    fn count_secrets() -> usize { 512 }
    fn count_shares() -> usize { 2187 }
    fn omega_secrets() -> u32 { 11641558 }
    fn omega_shares() -> u32 { 35260261 }
}

pub fn share_fft_sampling_sharing<C: Config, F>(b: &mut Bencher) 
where F: PrimeField + New<u32> + Encode<u32>, F::P: From<u32>, F::E: Clone
{
    let ref field = F::new(C::prime());
    let ref omega_secrets = field.encode(C::omega_secrets());
    let ref omega_shares = field.encode(C::omega_shares());
    
    let ref values = field.encode_slice(vec![5; C::count_secrets()]);
    
    b.iter(|| {
        let mut data = values.clone();
        ::numtheory::fft::fft2_inverse(field, &mut *data, omega_secrets);
        
        data.extend(vec![field.zero(); C::count_shares() - C::count_secrets()]);
        ::numtheory::fft::fft3(field, &mut *data, omega_shares);
    });
}

pub fn share_fft_sampling<C: Config, F>(b: &mut Bencher) 
where F: PrimeField + New<u32> + Encode<u32>, F::P: From<u32>, F::E: Clone
{
    let ref field = F::new(C::prime());
    let ref omega_secrets = field.encode(C::omega_secrets());
    
    let ref values = field.encode_slice(vec![5; C::count_secrets()]);
    
    b.iter(|| {
        let mut data = values.clone();
        ::numtheory::fft::fft2_inverse(field, &mut *data, &omega_secrets);
        
        let _shares = (1..C::count_shares() as u32+1)
            .map(|p| ::numtheory::mod_evaluate_polynomial(&data, &field.encode(p), field))
            .collect::<Vec<_>>();
    });
}

pub fn share_newton<C: Config, F>(b: &mut Bencher) 
where F: PrimeField + New<u32> + Encode<u32>, F::P: From<u32>, F::E: Clone
{
    let ref field = F::new(C::prime());
    
    let ref values = field.encode_slice(vec![5; C::count_secrets()]);
    let ref points = (1..C::count_secrets() as u32+1)
        .map(|p| field.sub(field.zero(), field.encode(p)))
        .collect::<Vec<_>>();
    
    b.iter(|| {
        let _shares = (1..C::count_shares() as u32+1)
            .map(|p| ::numtheory::newton_interpolation_at_point(&field.encode(p), points, values, field))
            .collect::<Vec<_>>();
    });
}

pub fn share_newton_pre<C: Config, F>(b: &mut Bencher) 
where F: PrimeField + New<u32> + Encode<u32>, F::P: From<u32>, F::E: Clone
{
    let ref field = F::new(C::prime());
    
    let ref values = field.encode_slice(vec![5; C::count_secrets()]);
    let ref points = (1..C::count_secrets() as u32+1)
        .map(|p| field.sub(field.zero(), field.encode(p)))
        .collect::<Vec<_>>();
    
    b.iter(|| {
        // this cannot be precomputed since it depends on the values
        let poly = ::numtheory::NewtonPolynomial::compute(points, values, field);
        
        let _shares = (1..C::count_shares() as u32+1)
            .map(|p| poly.evaluate(field.encode(p), field))
            .collect::<Vec<_>>();
    });
}

pub fn share_lagrange<C: Config, F>(b: &mut Bencher) 
where F: PrimeField + New<u32> + Encode<u32>, F::P: From<u32>, F::E: Clone
{
    let ref field = F::new(C::prime());
    
    let ref values = field.encode_slice(vec![5; C::count_secrets()]);
    let ref points = (1..C::count_secrets() as u32+1)
        .map(|p| field.sub(field.zero(), field.encode(p)))
        .collect::<Vec<_>>();
    
    b.iter(|| {
        let _shares = (1..C::count_shares() as u32+1)
            .map(|p| ::numtheory::lagrange_interpolation_at_point(&field.encode(p), points, values, field))
            .collect::<Vec<_>>();
    });
}

pub fn share_lagrange_pre<C: Config, F>(b: &mut Bencher) 
where F: PrimeField + New<u32> + Encode<u32>, F::P: From<u32>, F::E: Clone
{
    let ref field = F::new(C::prime());
    
    let ref values = field.encode_slice(vec![5; C::count_secrets()]);
    let ref points = (1..C::count_secrets() as u32+1)
        .map(|p| field.sub(field.zero(), field.encode(p)))
        .collect::<Vec<_>>();
    
    // this could be precomputed since it's independent of the values and the points are fixed
    let ref constants = (1..C::count_shares() as u32+1)
        .map(|p| ::numtheory::LagrangeConstants::compute(&field.encode(p), points, field))
        .collect::<Vec<_>>();
    
    b.iter(|| {
        let _shares = constants.iter()
            .map(|c| c.interpolate(values, field))
            .collect::<Vec<_>>();
    });
}

benchmark_group!(small
    , share_fft_sampling_sharing <Small, MontgomeryField32>
    , share_fft_sampling         <Small, MontgomeryField32>
    , share_newton               <Small, MontgomeryField32>
    , share_newton_pre           <Small, MontgomeryField32>
    , share_lagrange             <Small, MontgomeryField32>
    , share_lagrange_pre         <Small, MontgomeryField32>
);
 
benchmark_group!(medium
    , share_fft_sampling_sharing <Medium, MontgomeryField32>
    , share_fft_sampling         <Medium, MontgomeryField32>
    // , share_newton               <Medium, MontgomeryField32>
    , share_newton_pre           <Medium, MontgomeryField32>
    // , share_lagrange             <Medium, MontgomeryField32>
    , share_lagrange_pre         <Medium, MontgomeryField32>
);

benchmark_group!(large
    , share_fft_sampling_sharing <Large, MontgomeryField32>
    , share_lagrange_pre         <Large, MontgomeryField32>
);

benchmark_main!(
    // small
    // medium
    large
    // small, medium
);
