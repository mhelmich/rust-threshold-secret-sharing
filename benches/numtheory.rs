// Copyright (c) 2017 rust-threshold-secret-sharing developers

#[macro_use]
extern crate bencher;
extern crate threshold_secret_sharing as tss;

use bencher::Bencher;
use tss::*;

pub fn bench_binary_egcd(b: &mut Bencher) {
    b.iter(|| {
        let _ = tss::numtheory::binary_egcd(18461641736, 171);
    })
}

pub fn bench_euclidean_egcd(b: &mut Bencher) {
    b.iter(|| {
        let _ = tss::numtheory::gcd(18461641736, 171);
    })
}

benchmark_group!(egcd,
                 bench_euclidean_egcd,
                 bench_binary_egcd);


pub fn bench_weighted_sum_two_step(b: &mut Bencher) {
    let ref field = MontgomeryField32::new(746497_u32.into());
    let ref values = field.encode_slice(vec![5; 100]);
    let ref weights = field.encode_slice(vec![2; 100]);
    
    b.iter(|| {
        let _ = values.iter()
            .zip(weights)
            .map(|(v, w)| field.mul(v, w))
            .fold(field.zero(), |sum, term| field.add(sum, term));
    })
}

pub fn bench_weighted_sum_one_step(b: &mut Bencher) {
    let ref field = MontgomeryField32::new(746497_u32.into());
    let ref values = field.encode_slice(vec![5; 100]);
    let ref weights = field.encode_slice(vec![2; 100]);
    
    b.iter(|| {
        let _ = values.iter()
            .zip(weights)
            .fold(field.zero(), |sum, (v, w)| field.add(sum, field.mul(v, w)));
    })
}

pub fn bench_weighted_sum_for(b: &mut Bencher) {
    let ref field = MontgomeryField32::new(746497_u32.into());
    let ref values = field.encode_slice(vec![5; 100]);
    let ref weights = field.encode_slice(vec![2; 100]);
    
    b.iter(|| {
        let mut sum = field.zero();
        for i in 0..values.len() {
            let v = values[i];
            let w = weights[i];
            sum = field.add(sum, field.mul(v, w));
        }
    })
}

benchmark_group!(weighted_sum,
                 bench_weighted_sum_two_step,
                 bench_weighted_sum_one_step,
                 bench_weighted_sum_for);

benchmark_main!(egcd, weighted_sum);
