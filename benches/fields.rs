// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.

#[macro_use]
extern crate bencher;
extern crate threshold_secret_sharing as tss;

use bencher::Bencher;
use tss::*;

pub fn bench_add<F>(bencher: &mut Bencher)
where
    F: PrimeField + New<u32> + Encode<u32>,
    F::P: From<u32>,
{
    let field = F::new(746497);

    let a = field.encode(1231231);
    let b = field.encode(423421);

    bencher.iter(|| {
        let _ = field.add(&a, &b);
    });
}

pub fn bench_mul<F>(bencher: &mut Bencher)
where
    F: PrimeField + New<u32> + Encode<u32>,
    F::P: From<u32>,
{
    let field = F::new(746497);

    let a = field.encode(1231231);
    let b = field.encode(423421);

    bencher.iter(|| {
        let _ = field.mul(&a, &b);
    });
}

benchmark_group!(
    add,
    bench_add<NaturalPrimeField<i64>>,
    bench_add<MontgomeryField32>,
    bench_add<LargePrimeField>
);

benchmark_group!(
    mul,
    bench_mul<NaturalPrimeField<i64>>,
    bench_mul<MontgomeryField32>,
    bench_mul<LargePrimeField>
);

benchmark_main!(add, mul);
