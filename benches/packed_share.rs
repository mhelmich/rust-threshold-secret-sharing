// Copyright (c) 2017 rust-threshold-secret-sharing developers

#![allow(dead_code)]

#[macro_use]
extern crate bencher;
extern crate threshold_secret_sharing as tss;

use bencher::Bencher;
use tss::*;

pub trait Config {
    fn privacy() -> usize;
    fn secrets() -> usize;
    fn shares() -> usize;
    fn prime() -> u32;
    fn order_secrets() -> usize;
    fn omega_secrets() -> u32;
    fn order_shares() -> usize;
    fn omega_shares() -> u32;
}

struct Small;
impl Config for Small {
    fn privacy() -> usize {
        2
    }
    fn secrets() -> usize {
        4
    } // 2 + 4 + 1 <= 8
    fn shares() -> usize {
        8
    }
    fn prime() -> u32 {
        433
    }
    fn order_secrets() -> usize {
        8
    }
    fn omega_secrets() -> u32 {
        354
    }
    fn order_shares() -> usize {
        9
    }
    fn omega_shares() -> u32 {
        150
    }
}

struct Medium;
impl Config for Medium {
    fn privacy() -> usize {
        7
    }
    fn secrets() -> usize {
        8
    } // 7 + 8 + 1 <= 16
    fn shares() -> usize {
        26
    }
    fn prime() -> u32 {
        746497
    }
    fn order_secrets() -> usize {
        16
    }
    fn omega_secrets() -> u32 {
        328238
    }
    fn order_shares() -> usize {
        27
    }
    fn omega_shares() -> u32 {
        514357
    }
}

struct Large;
impl Config for Large {
    fn privacy() -> usize {
        20
    }
    fn secrets() -> usize {
        40
    } // 20 + 40 + 1 <= 64
    fn shares() -> usize {
        80
    }
    fn prime() -> u32 {
        746497
    }
    fn order_secrets() -> usize {
        64
    }
    fn omega_secrets() -> u32 {
        181622
    }
    fn order_shares() -> usize {
        81
    }
    fn omega_shares() -> u32 {
        69177
    }
}

struct Huge;
impl Config for Huge {
    fn privacy() -> usize {
        61
    }
    fn secrets() -> usize {
        66
    } // 61 + 66 + 1 <= 128
    fn shares() -> usize {
        242
    }
    fn prime() -> u32 {
        746497
    }
    fn order_secrets() -> usize {
        128
    }
    fn omega_secrets() -> u32 {
        275374
    }
    fn order_shares() -> usize {
        243
    }
    fn omega_shares() -> u32 {
        595577
    }
}

pub fn share_fft_sampling_sharing<C: Config, F>(b: &mut Bencher)
where
    F: PrimeField + New<u32> + Encode<u32>,
    F::P: From<u32>,
    F::E: Clone,
{
    let ref field = F::new(C::prime());
    let ref omega_secrets = field.encode(C::omega_secrets());
    let ref omega_shares = field.encode(C::omega_shares());

    let ref values = field.encode_slice(vec![5; C::order_secrets()]);

    b.iter(|| {
        let mut data = values.clone();
        ::numtheory::fft::fft2_inverse(field, &mut *data, omega_secrets);

        data.extend(vec![field.zero(); C::order_shares() - C::order_secrets()]);
        ::numtheory::fft::fft3(field, &mut *data, omega_shares);
    });
}

pub fn share_fft_sampling<C: Config, F>(b: &mut Bencher)
where
    F: PrimeField + New<u32> + Encode<u32>,
    F::P: From<u32>,
    F::E: Clone,
{
    let ref field = F::new(C::prime());
    let ref omega_secrets = field.encode(C::omega_secrets());

    let ref values = field.encode_slice(vec![5; C::order_secrets()]);

    b.iter(|| {
        let mut data = values.clone();
        ::numtheory::fft::fft2_inverse(field, &mut *data, &omega_secrets);

        let _shares = (1..C::shares() as u32 + 1)
            .map(|p| ::numtheory::mod_evaluate_polynomial(&data, &field.encode(p), field))
            .collect::<Vec<_>>();
    });
}

pub fn share_newton<C: Config, F>(b: &mut Bencher)
where
    F: PrimeField + New<u32> + Encode<u32>,
    F::P: From<u32>,
    F::E: Clone,
{
    let ref field = F::new(C::prime());

    let ref values = field.encode_slice(vec![5; C::privacy() + C::secrets()]);
    let ref points = (1..C::privacy() + C::secrets() + 1)
        .map(|p| field.sub(field.zero(), field.encode(p as u32)))
        .collect::<Vec<_>>();

    b.iter(|| {
        let _shares = (1..C::shares() as u32 + 1)
            .map(|p| {
                ::numtheory::newton_interpolation_at_point(&field.encode(p), points, values, field)
            })
            .collect::<Vec<_>>();
    });
}

pub fn share_newton_pre<C: Config, F>(b: &mut Bencher)
where
    F: PrimeField + New<u32> + Encode<u32>,
    F::P: From<u32>,
    F::E: Clone,
{
    let ref field = F::new(C::prime());

    let ref values = field.encode_slice(vec![5; C::privacy() + C::secrets()]);
    let ref points = (1..C::privacy() + C::secrets() + 1)
        .map(|p| field.sub(field.zero(), field.encode(p as u32)))
        .collect::<Vec<_>>();

    b.iter(|| {
        // this cannot be precomputed since it depends on the values
        let poly = ::numtheory::NewtonPolynomial::compute(points, values, field);

        let _shares = (1..C::shares() as u32 + 1)
            .map(|p| poly.evaluate(field.encode(p), field))
            .collect::<Vec<_>>();
    });
}

pub fn share_lagrange<C: Config, F>(b: &mut Bencher)
where
    F: PrimeField + New<u32> + Encode<u32>,
    F::P: From<u32>,
    F::E: Clone,
{
    let ref field = F::new(C::prime());

    let ref values = field.encode_slice(vec![5; C::privacy() + C::secrets()]);
    let ref points = (1..C::privacy() + C::secrets() + 1)
        .map(|p| field.sub(field.zero(), field.encode(p as u32)))
        .collect::<Vec<_>>();

    b.iter(|| {
        let _shares = (1..C::shares() as u32 + 1)
            .map(|p| {
                ::numtheory::lagrange_interpolation_at_point(
                    &field.encode(p),
                    points,
                    values,
                    field,
                )
            })
            .collect::<Vec<_>>();
    });
}

pub fn share_lagrange_pre<C: Config, F>(b: &mut Bencher)
where
    F: PrimeField + New<u32> + Encode<u32>,
    F::P: From<u32>,
    F::E: Clone,
{
    let ref field = F::new(C::prime());

    let ref values = field.encode_slice(vec![5; C::privacy() + C::secrets()]);
    let ref points = (1..C::privacy() + C::secrets() + 1)
        .map(|p| field.sub(field.zero(), field.encode(p as u32)))
        .collect::<Vec<_>>();

    // this could be precomputed since it's independent of the values and the points are fixed
    let ref constants = (1..C::shares() as u32 + 1)
        .map(|p| ::numtheory::LagrangeConstants::compute(&field.encode(p), points, field))
        .collect::<Vec<_>>();

    b.iter(|| {
        let _shares = constants
            .iter()
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
    , share_newton               <Medium, MontgomeryField32>
    , share_newton_pre           <Medium, MontgomeryField32>
    , share_lagrange             <Medium, MontgomeryField32>
    , share_lagrange_pre         <Medium, MontgomeryField32>
);

benchmark_group!(large
    , share_fft_sampling_sharing <Large, MontgomeryField32>
    , share_fft_sampling         <Large, MontgomeryField32>
    , share_newton_pre           <Large, MontgomeryField32>
    , share_lagrange_pre         <Large, MontgomeryField32>
);

benchmark_group!(huge
    , share_fft_sampling_sharing <Huge, MontgomeryField32>
    , share_fft_sampling         <Huge, MontgomeryField32>
    , share_newton_pre           <Huge, MontgomeryField32>
    , share_lagrange_pre         <Huge, MontgomeryField32>
);

benchmark_main!(small, medium, large, huge);
