fn main() {
    let mut link_type = "dylib";
    println!("cargo:rerun-if-env-changed=GF2POLY_LIB_PATH");
    println!("cargo:rerun-if-env-changed=GF2POLY_STATIC_LIB");
    for (key, value) in std::env::vars() {
        match key.as_str() {
            "GF2POLY_STATIC_LIB" => {
                if value == "1" {
                    link_type = "static";
                }
            }
            "GF2POLY_LIB_PATH" => {
                println!("cargo:rustc-link-search=native={}", value);
            }
            _ => continue,
        }
    }
    println!("cargo:rustc-link-lib={link_type}=gf2x");
}
