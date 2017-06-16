// Copyright (c) 2017 rust-threshold-secret-sharing developers

extern crate threshold_secret_sharing as tss;
extern crate rand;

#[cfg(not(
    all(feature="safety_override",
        any(feature="largefield_ramp", feature="largefield_gmp")
    )))]
fn main() {
    println!("Please run with '--features \"largefield safety_override\"'");
}

#[cfg(
    all(feature="safety_override",
        any(feature="largefield_ramp", feature="largefield_gmp")
    ))]
#[allow(non_snake_case)]
fn main() {
    
    use tss::*;

    let n = 1024;
    let m = 2187;
    
    let field = LargePrimeField::new("2168493841578655774908481580141050902529");
    let omega_n = field.encode("575568907032575917226189174489221138041");
    let omega_m = field.encode("611000397685398983853362438068827993522");
    
    let n_A = n/2;
    let n_B = n;
    let omega_secrets_A = field.pow(&omega_n, 2);
    let omega_secrets_B = omega_n.clone();
    let omega_shares = omega_m.clone();
    
    let secret_count = 250;
    
    // sanity checks
    use std::collections::HashSet;
    // - roots of unity
    assert![ field.pow(&omega_secrets_A, n_A) == field.one() ];
    assert![ field.pow(&omega_secrets_B, n_B) == field.one() ];
    assert![ field.pow(&omega_shares, m) == field.one() ];
    // - principal
    let group_n_A: HashSet<_> = (1..n_A).map(|e| field.pow(&omega_secrets_A, e) ).collect();
    let group_n_B: HashSet<_> = (1..n_B).map(|e| field.pow(&omega_secrets_B, e) ).collect();
    let group_m: HashSet<_> = (1..m).map(|e| field.pow(&omega_shares, e) ).collect();
    assert![ group_n_A.len() as u32 == n_A - 1 ];
    assert![ group_n_B.len() as u32 == n_B - 1 ];
    assert![ group_m.len() as u32 == m - 1 ];
    assert![ ! group_n_A.contains(&field.one()) ];
    assert![ ! group_n_B.contains(&field.one()) ];
    assert![ ! group_m.contains(&field.one()) ];
    // - points overlap
    assert![ group_n_A.intersection(&group_m).count() == 0 ];
    assert![ group_n_B.intersection(&group_m).count() == 0 ];
    assert![ group_n_A.intersection(&group_n_B).count() as u32 == n/2 - 1 ];
    
    let pss_A = PackedSecretSharing {
        secret_count: secret_count,
        threshold: n_A as usize - secret_count - 1,
        share_count: m as usize - 1,
        omega_secrets: omega_secrets_A,
        omega_shares: omega_shares.clone(),
        field: field.clone(),
    };
    
    let pss_B = PackedSecretSharing {
        secret_count: secret_count,
        threshold: n_B as usize - secret_count - 1,
        share_count: m as usize - 1,
        omega_secrets: omega_secrets_B,
        omega_shares: omega_shares.clone(),
        field: field.clone(),
    };
    
    println!("Sharing for A..");
    let secrets_A = field.encode_slice((0..secret_count as u32).collect::<Vec<_>>());
    println!("{:?}", secrets_A);
    let shares_A = pss_A.share(&secrets_A);
    
    println!("Sharing for B..");
    let secrets_B = field.encode_slice((0..secret_count as u32).collect::<Vec<_>>());
    println!("{:?}", secrets_B);
    let randomness_B = {
        use rand;
        let mut rng = rand::OsRng::new().unwrap();
        field.sample_with_replacement(pss_B.threshold, &mut rng)
    };
    let mut secrets_and_randomess_B = randomness_B;
    for (i,s) in secrets_B.into_iter().enumerate() {
        secrets_and_randomess_B.insert(2*i+1, s);
    }
    let shares_B = pss_B.deterministic_share(&secrets_and_randomess_B);
    
    println!("Summing..");
    let shares_sum = shares_A.iter().zip(&shares_B).map(|(a, b)| field.add(a, b)).collect::<Vec<_>>();
    
    println!("Reconstructing..");
    let indices = (0..shares_sum.len() as u32).collect::<Vec<_>>();
    let secrets_sum = pss_B.fully_reconstruct(&*indices, &shares_sum).into_iter()
        .enumerate()
        .filter(|&(i,_)| i % 2 == 1)
        .map(|(_,s)| s)
        .take(secret_count)
        .collect::<Vec<_>>();
    println!("{:?}", secrets_sum);
}
