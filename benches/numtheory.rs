// Copyright (c) 2017 rust-threshold-secret-sharing developers

#[macro_use]
extern crate bencher;
extern crate threshold_secret_sharing as tss;

use bencher::Bencher;

pub fn bench_binary_egcd(b: &mut Bencher) {
    use tss::numtheory::binary_egcd;
    b.iter(|| {
        let _ = binary_egcd(18461641736, 171);
    })
}

pub fn bench_euclidean_egcd(b: &mut Bencher) {
    use tss::numtheory::gcd;
    b.iter(|| {
        let _ = gcd(18461641736, 171);
    })
}

benchmark_group!(egcd,
                 bench_euclidean_egcd,
                 bench_binary_egcd);

benchmark_main!(egcd);
