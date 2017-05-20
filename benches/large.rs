// Copyright (c) 2017 rust-threshold-secret-sharing developers

#[macro_use]
extern crate bencher;
extern crate threshold_secret_sharing as tss;
extern crate ramp;

use tss::*;

use bencher::Bencher;
use tss::packed::*;

trait Size {
    fn prime() -> ramp::Int;
}

pub fn bench_100_packed<P: Prime>(b: &mut Bencher) 
{
    let prime = P::prime();
    let field = LargePrimeField::new(prime);
    
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

benchmark_group!(shamir_vs_packed,
     shamir_vs_packed::bench_100_shamir<MontgomeryField32>,
     shamir_vs_packed::bench_100_packed<MontgomeryField32>,
     shamir_vs_packed::bench_100_shamir<NaturalPrimeField<i64>>,
     shamir_vs_packed::bench_100_packed<NaturalPrimeField<i64>>,
     shamir_vs_packed::bench_100_shamir<LargePrimeField>,
     shamir_vs_packed::bench_100_packed<LargePrimeField>
 );

benchmark_main!(shamir_vs_packed, packed);
