// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.

//! This module implements in-place 2-radix and 3-radix numeric theory
//! transformations (FFT on modular fields) by in-place Cooley-Tukey algorithms.

use fields::Field;
use fields::Encode;

/// 2-radix FFT.
///
/// * zp is the modular field
/// * data is the data to transform
/// * omega is the root-of-unity to use
///
/// `data.len()` must be a power of 2. omega must be a root of unity of order
/// `data.len()`
pub fn fft2<F>(zp: &F, data: &mut [F::E], omega: F::E) 
where F: Field, F::E: Clone
{
    fft2_in_place_rearrange(zp, &mut *data);
    fft2_in_place_compute(zp, &mut *data, omega);
}

/// 2-radix inverse FFT.
///
/// * zp is the modular field
/// * data is the data to transform
/// * omega is the root-of-unity to use
///
/// `data.len()` must be a power of 2. omega must be a root of unity of order
/// `data.len()`
pub fn fft2_inverse<F>(zp: &F, data: &mut [F::E], omega: F::E) 
where F: Field + Encode<u32>, F::E: Clone
{
    let omega_inv = zp.inv(omega);
    let len = data.len();
    let len_inv = zp.inv(zp.encode(len as u32));
    fft2(zp, data, omega_inv);
    for mut x in data {
        *x = zp.mul(&*x, &len_inv);
    }
}

fn fft2_in_place_rearrange<F>(_zp: &F, data: &mut [F::E]) 
where F: Field
{
    let mut target = 0;
    for pos in 0..data.len() {
        if target > pos {
            data.swap(target, pos)
        }
        let mut mask = data.len() >> 1;
        while target & mask != 0 {
            target &= !mask;
            mask >>= 1;
        }
        target |= mask;
    }
}

fn fft2_in_place_compute<F>(zp: &F, data: &mut [F::E], omega: F::E) 
where F: Field, F::E: Clone
{
    let mut depth = 0usize;
    while 1usize << depth < data.len() {
        let step = 1usize << depth;
        let jump = 2 * step;
        let factor_stride = zp.pow(&omega, (data.len() / step / 2) as u32);
        let mut factor = zp.one();
        for group in 0usize..step {
            let mut pair = group;
            while pair < data.len() {
                let x = data[pair].clone();
                let y = zp.mul(&data[pair + step], &factor);

                data[pair] = zp.add(&x, &y);
                data[pair + step] = zp.sub(&x, &y);

                pair += jump;
            }
            factor = zp.mul(factor, &factor_stride);
        }
        depth += 1;
    }
}

/// 3-radix FFT.
///
/// * zp is the modular field
/// * data is the data to transform
/// * omega is the root-of-unity to use
///
/// `data.len()` must be a power of 2. omega must be a root of unity of order
/// `data.len()`
pub fn fft3<F>(zp: &F, data: &mut [F::E], omega: F::E) 
where F: Field, F::E: Clone
{
    fft3_in_place_rearrange(zp, &mut *data);
    fft3_in_place_compute(zp, &mut *data, omega);
}

/// 3-radix inverse FFT.
///
/// * zp is the modular field
/// * data is the data to transform
/// * omega is the root-of-unity to use
///
/// `data.len()` must be a power of 2. omega must be a root of unity of order
/// `data.len()`
pub fn fft3_inverse<F>(zp: &F, data: &mut [F::E], omega: F::E) 
where F: Field + Encode<u32>, F::E: Clone
{
    let omega_inv = zp.inv(omega);
    let len_inv = zp.inv(zp.encode(data.len() as u32));
    fft3(zp, data, omega_inv);
    for mut x in data {
        *x = zp.mul(&*x, &len_inv);
    }
}

fn trigits_len(n: usize) -> usize {
    let mut result = 1;
    let mut value = 3;
    while value < n + 1 {
        result += 1;
        value *= 3;
    }
    result
}

fn fft3_in_place_rearrange<F: Field>(_zp: &F, data: &mut [F::E]) {
    let mut target = 0isize;
    let trigits_len = trigits_len(data.len() - 1);
    let mut trigits: Vec<u8> = ::std::iter::repeat(0).take(trigits_len).collect();
    let powers: Vec<isize> = (0..trigits_len).map(|x| 3isize.pow(x as u32)).rev().collect();
    for pos in 0..data.len() {
        if target as usize > pos {
            data.swap(target as usize, pos)
        }
        for pow in 0..trigits_len {
            if trigits[pow] < 2 {
                trigits[pow] += 1;
                target += powers[pow];
                break;
            } else {
                trigits[pow] = 0;
                target -= 2 * powers[pow];
            }
        }
    }
}

fn fft3_in_place_compute<F>(zp: &F, data: &mut [F::E], omega: F::E) 
where F: Field, F::E: Clone
{
    let mut step = 1;
    let big_omega = zp.pow(&omega, (data.len() as u32 / 3));
    let big_omega_sq = zp.mul(&big_omega, &big_omega);
    while step < data.len() {
        let jump = 3 * step;
        let factor_stride = zp.pow(&omega, (data.len() / step / 3) as u32);
        let mut factor = zp.one();
        for group in 0usize..step {
            let factor_sq = zp.mul(&factor, &factor);
            let mut pair = group;
            while pair < data.len() {
                let x = data[pair].clone();
                let y = zp.mul(&data[pair + step], &factor);
                let z = zp.mul(&data[pair + 2 * step], &factor_sq);

                data[pair] = zp.add(zp.add(&x, &y), &z);
                data[pair + step] = zp.add(zp.add(&x, zp.mul(&big_omega, &y)), zp.mul(&big_omega_sq, &z));
                data[pair + 2 * step] = zp.add(zp.add(&x, zp.mul(&big_omega_sq, &y)), zp.mul(&big_omega, &z));

                pair += jump;
            }
            factor = zp.mul(&factor, &factor_stride);
        }
        step = jump;
    }
}

#[cfg(test)]
pub mod test {
    use super::*;
    use ::fields::{Field, PrimeField, Encode, Decode};
    
    pub fn from<F>(field: &F, data: &[u32]) -> Vec<F::E> 
    where F: Field + Encode<u32>, F::E: Clone
    {
        data.iter().map(|&x| field.encode(x)).collect()
    }

    pub fn back<F>(field: &F, data: &[F::E]) -> Vec<u32> 
    where F: Field + Decode<u32>, F::E: Clone
    {
        data.iter().map(|x| field.decode(x)).collect()
    }

    pub fn test_fft2<F>()
    where F: PrimeField + Encode<u32> + Decode<u32>, F::E: Clone, F::P: From<u32>
    {
        // field is Z_433 in which 354 is an 8th root of unity
        let zp = F::new(433.into());
        let omega = zp.encode(354);

        let mut data = from(&zp, &[1, 2, 3, 4, 5, 6, 7, 8]);
        fft2(&zp, &mut data, omega);
        assert_eq!(back(&zp, &data), [36, 303, 146, 3, 429, 422, 279, 122]);
    }

    pub fn test_fft2_inverse<F>() 
    where F: PrimeField + Encode<u32> + Decode<u32>, F::E: Clone, F::P: From<u32>
    {
        // field is Z_433 in which 354 is an 8th root of unity
        let zp = F::new(433.into());
        let omega = zp.encode(354);

        let mut data = from(&zp, &[36, 303, 146, 3, 429, 422, 279, 122]);
        fft2_inverse(&zp, &mut *data, omega);
        assert_eq!(back(&zp, &data), [1, 2, 3, 4, 5, 6, 7, 8])
    }

    pub fn test_fft2_big<F>() 
    where F: PrimeField + Encode<u32> + Decode<u32>, F::E: Clone, F::P: From<u32>
    {
        let zp = F::new(5038849.into());
        let omega = zp.encode(4318906);

        let mut data: Vec<_> = (0..256).map(|a| zp.encode(a)).collect();
        fft2(&zp, &mut *data, omega.clone());
        fft2_inverse(&zp, &mut data, omega.clone());

        assert_eq!(back(&zp, &data), (0..256).collect::<Vec<_>>());
    }

    pub fn test_fft3<F>() 
    where F: PrimeField + Encode<u32> + Decode<u32>, F::E: Clone, F::P: From<u32>
    {
        // field is Z_433 in which 150 is an 9th root of unity
        let zp = F::new(433.into());
        let omega = zp.encode(150);

        let mut data = from(&zp, &[1, 2, 3, 4, 5, 6, 7, 8, 9]);
        fft3(&zp, &mut data, omega);
        assert_eq!(back(&zp, &data), [45, 404, 407, 266, 377, 47, 158, 17, 20]);
    }

    pub fn test_fft3_inverse<F>() 
    where F: PrimeField + Encode<u32> + Decode<u32>, F::E: Clone, F::P: From<u32>
    {
        // field is Z_433 in which 150 is an 9th root of unity
        let zp = F::new(433.into());
        let omega = zp.encode(150);

        let mut data = from(&zp, &[45, 404, 407, 266, 377, 47, 158, 17, 20]);
        fft3_inverse(&zp, &mut *data, omega);
        assert_eq!(back(&zp, &data), [1, 2, 3, 4, 5, 6, 7, 8, 9])
    }

    pub fn test_fft3_big<F>() 
    where F: PrimeField + Encode<u32> + Decode<u32>, F::E: Clone, F::P: From<u32>
    {
        let zp = F::new(5038849.into());
        let omega = zp.encode(1814687);

        let mut data: Vec<_> = (0..19683).map(|a| zp.encode(a)).collect();
        fft3(&zp, &mut data, omega.clone());
        fft3_inverse(&zp, &mut data, omega.clone());

        assert_eq!(back(&zp, &data), (0..19683).collect::<Vec<_>>());
    }
}
