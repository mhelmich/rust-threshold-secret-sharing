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

use tss::*;

mod shamir_vs_packed {

    use super::*;
    use bencher::Bencher;

    pub fn bench_100_shamir<F>(b: &mut Bencher) 
    where F: PrimeField + Clone + Encode<u32>, F::P: From<u32>, F::E: Clone
    {
        let field = F::new(746497.into());
        
        let tss = ShamirSecretSharing {
            threshold: 155 / 3,
            share_count: 728 / 3,
            field: field.clone(),
        };

        let all_secrets = field.encode_slice(vec![5 ; 100]);
            
        b.iter(|| {
            let _shares: Vec<Vec<_>> = all_secrets.iter()
                .map(|secret| tss.share(secret.clone()))
                .collect();
        });
    }

    pub fn bench_100_packed<F>(b: &mut Bencher) 
    where F: PrimeField + Clone + Encode<u32>, F::P: From<u32>, F::E: Clone
    {
        let field = F::new(746497.into());
        
        let pss = PackedSecretSharing {
            threshold: 155,
            share_count: 728,
            secret_count: 100,
            omega_secrets: field.encode(95660),
            omega_shares: field.encode(610121),
            field: field.clone(),
        };
        
        let all_secrets = field.encode_slice(vec![5 ; 100]);
        
        b.iter(|| {
            let _shares = pss.share(&all_secrets);
        })
    }

}

#[cfg(not(feature="largefield"))]
benchmark_group!(shamir_vs_packed,
                 shamir_vs_packed::bench_100_shamir<MontgomeryField32>,
                 shamir_vs_packed::bench_100_packed<MontgomeryField32>,
                 shamir_vs_packed::bench_100_shamir<NaturalPrimeField<i64>>,
                 shamir_vs_packed::bench_100_packed<NaturalPrimeField<i64>>
 );

#[cfg(feature="largefield")]
benchmark_group!(shamir_vs_packed,
                 shamir_vs_packed::bench_100_shamir<MontgomeryField32>,
                 shamir_vs_packed::bench_100_packed<MontgomeryField32>,
                 shamir_vs_packed::bench_100_shamir<NaturalPrimeField<i64>>,
                 shamir_vs_packed::bench_100_packed<NaturalPrimeField<i64>>,
                 shamir_vs_packed::bench_100_shamir<LargePrimeField>,
                 shamir_vs_packed::bench_100_packed<LargePrimeField>
 );

benchmark_main!(shamir_vs_packed);
