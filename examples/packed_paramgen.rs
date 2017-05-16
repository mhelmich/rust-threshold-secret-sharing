// Copyright (c) 2017 rust-threshold-secret-sharing developers

extern crate threshold_secret_sharing as tss;

#[cfg(not(feature="paramgen"))]
fn main() {}

#[cfg(feature="paramgen")]
fn main() {

    use tss::packed::PackedSecretSharing;

    let threshold = 5; // privacy threshold
    let secret_count = 10; // number of secrets packed together
    // the sum of the above two must be a power of two minus one: 5 + 10 = 2^4 - 1
    
    let share_count = 26; // total number of shares
    // this must be a power of three minus one: 26 = 3^3 - 1
    
    let min_size = 1 << 16;
    // we want a prime field supporting 16 bit numbers

    let ref pss = PackedSecretSharing::new_with_min_size(
        threshold,
        secret_count,
        share_count,
        min_size
    );
    println!("{:?}", pss);
    
    // let ref expected_pss = PackedSecretSharing { 
    //     threshold: 5,
    //     secret_count: 10,
    //     share_count: 26,
    //     prime: 66529,
    //     omega_secrets: 17871,
    //     omega_shares: 994
    // };
    // assert_eq!(pss, expected_pss);
    
    let secrets = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
    let shares = pss.share(&secrets);
    
    let indices: Vec<usize> = (0..shares.len()).collect();
    let reconstructed_secrets = pss.reconstruct(&indices, &shares);
    
    // assert_eq!(secrets, reconstructed_secrets);
    // above will likely fail since some values may be `x - prime` instead of `x`
    
    assert_eq!(secrets, tss::positivise(&reconstructed_secrets, pss.prime));
    // the above ensures a canonical representation of numbers in `[0, prime)`
    
}
