//! Compiles the Simulink-generated environment C code into this crate and
//! tells rustc where luna's source lives.
//!
//! Environment variables (all optional):
//! - `HOPPER_ENV_CODE`: folder with the generated code (`hopper_env.mk`,
//!   `hopper_env.c`, ...). Default: `../build/hopper_env_ert_rtw`.
//! - `MATLAB_ROOT`: MATLAB install, for the generated code's headers.
//!   Default: `C:/Program Files/MATLAB/R2026a`.
//! - `HOPPER_CC`: C compiler. Default: `gcc` on PATH. Use the same MinGW
//!   MATLAB builds with. `ar` is taken from the same folder.
//! - `LUNA_DIR`: luna checkout. Default: `luna` next to the hopper repo.

use std::{
    env, fs,
    path::{Path, PathBuf},
};

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());

    let code = env_path("HOPPER_ENV_CODE", manifest.join("../build/hopper_env_ert_rtw"));
    let matlab = env_path("MATLAB_ROOT", PathBuf::from("C:/Program Files/MATLAB/R2026a"));
    // sil -> fsw_bridge -> hopper_6dof_sizer -> sims -> hopper -> (GitHub folder)
    let luna = env_path("LUNA_DIR", manifest.join("../../../../../luna"));

    let makefile = code.join("hopper_env.mk");
    let mk = fs::read_to_string(&makefile).unwrap_or_else(|e| {
        panic!(
            "cannot read {} ({e}); run build_env_c in MATLAB first or set HOPPER_ENV_CODE",
            makefile.display()
        )
    });

    let mut build = cc::Build::new();
    if let Ok(compiler) = env::var("HOPPER_CC") {
        let compiler = PathBuf::from(compiler);
        if let Some(bin) = compiler.parent() {
            build.archiver(bin.join("ar.exe"));
        }
        build.compiler(compiler);
    }

    // Same sources and -D flags MATLAB's own build used
    for src in make_list(&mk, "SRCS") {
        let file = code.join(Path::new(&src).file_name().unwrap());
        println!("cargo:rerun-if-changed={}", file.display());
        build.file(file);
    }
    for line in mk.lines().filter(|l| l.starts_with("DEFINES_")) {
        for flag in line.split_whitespace().filter(|t| t.starts_with("-D")) {
            match flag[2..].split_once('=') {
                Some((name, value)) => build.define(name, value),
                None => build.define(&flag[2..], None),
            };
        }
    }

    build
        .file(manifest.join("csrc/env_shim.c"))
        .include(&code)
        .include(code.parent().unwrap())
        .include(matlab.join("extern/include"))
        .include(matlab.join("simulink/include"))
        .include(matlab.join("rtw/c/src"))
        .warnings(false)
        .compile("hopper_env");

    println!("cargo:rerun-if-changed={}", makefile.display());
    println!("cargo:rerun-if-changed=csrc/env_shim.c");
    for var in ["HOPPER_ENV_CODE", "MATLAB_ROOT", "HOPPER_CC", "LUNA_DIR"] {
        println!("cargo:rerun-if-env-changed={var}");
    }

    // include!() wants forward slashes and no \\?\ prefix
    let luna_str = luna.to_string_lossy().replace('\\', "/");
    for file in ["common/src/comm/ctv.rs", "flight2/src/control/lqr.rs"] {
        let path = luna.join(file);
        assert!(path.is_file(), "luna file not found: {} (set LUNA_DIR)", path.display());
        println!("cargo:rerun-if-changed={}", path.display());
    }
    println!("cargo:rustc-env=LUNA_DIR={luna_str}");
}

fn env_path(var: &str, default: PathBuf) -> PathBuf {
    env::var(var).map(PathBuf::from).unwrap_or(default)
}

/// Values of a `NAME = a b c` makefile variable.
fn make_list(mk: &str, name: &str) -> Vec<String> {
    mk.lines()
        .find_map(|l| l.strip_prefix(name)?.trim_start().strip_prefix('='))
        .map(|rest| rest.split_whitespace().map(|s| s.replace("$(START_DIR)/", "")).collect())
        .unwrap_or_default()
}
