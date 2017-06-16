// Copyright (c) 2017 rust-threshold-secret-sharing developers

#[macro_use] extern crate bencher;
extern crate threshold_secret_sharing as tss;

use tss::*;
use bencher::Bencher;

pub trait Config {
    fn prime() -> &'static str;
    fn n() -> usize;
    fn m() -> usize;
    fn secret_count() -> usize;
    fn omega_n() -> &'static str;
    fn omega_m() -> &'static str;
}

struct Small;
impl Config for Small {
    fn prime() -> &'static str { "27384083395425770033880312552185387014657" }
    fn n() -> usize { 256 }
    fn m() -> usize { 729 }
    fn secret_count() -> usize { 250 }
    fn omega_n() -> &'static str { "11760837463486371896102955950879363236841" }
    fn omega_m() -> &'static str { "282515243551697216317750839613713827799" }
}

struct Medium;
impl Config for Medium {
    fn prime() -> &'static str { "11957119570273711187856305448866671914241" }
    fn n() -> usize { 256 }
    fn m() -> usize { 2187 }
    fn secret_count() -> usize { 250 }
    fn omega_n() -> &'static str { "3069700116802786132417224549683499080751" }
    fn omega_m() -> &'static str { "4240725152029417998256133660574413855628" }
}

struct Large;
impl Config for Large {
    fn prime() -> &'static str { "2168493841578655774908481580141050902529" }
    fn n() -> usize { 1024 }
    fn m() -> usize { 2187 }
    fn secret_count() -> usize { 250 }
    fn omega_n() -> &'static str { "575568907032575917226189174489221138041" }
    fn omega_m() -> &'static str { "611000397685398983853362438068827993522" }
}

pub fn bench_share<C: Config>(b: &mut Bencher) 
{
    let field = LargePrimeField::new(C::prime());
    let omega_n = field.encode(C::omega_n());
    let omega_m = field.encode(C::omega_m());
    let n = C::n();
    let m = C::m();
    let secret_count = C::secret_count();
    
    let pss = PackedSecretSharing {
        secret_count: secret_count,
        threshold: n - secret_count - 1,
        share_count: m - 1,
        omega_secrets: omega_n,
        omega_shares: omega_m,
        field: field.clone(),
    };
    
    let secrets = field.encode_slice((0..secret_count as u32).collect::<Vec<_>>());
    
    b.iter(|| {
        let _shares = pss.share(&secrets);
    })
    
    // 
    // println!("Summing..");
    // let shares_sum = shares_A.iter().zip(&shares_B).map(|(a, b)| field.add(a, b)).collect::<Vec<_>>();
    // 
    // println!("Reconstructing..");
    // let indices = (0..shares_sum.len() as u32).collect::<Vec<_>>();
    // let secrets_sum = pss_B.fully_reconstruct(&*indices, &shares_sum).into_iter()
    //     .enumerate()
    //     .filter(|&(i,_)| i % 2 == 1)
    //     .map(|(_,s)| s)
    //     .take(secret_count)
    //     .collect::<Vec<_>>();
    // println!("{:?}", secrets_sum);
    // 
}

benchmark_group!(group,
    bench_share<Small>,
    bench_share<Medium>,
    bench_share<Large>
);

benchmark_main!(group);
