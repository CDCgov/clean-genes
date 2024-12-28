use clap::Parser;
use std::path::Path;

#[derive(Parser, Default, Debug)]
#[command(
    name = "clean-genes",
    author = "David E. Hufnagel",
    version = env!("CARGO_PKG_VERSION"),
    about = "A CLI tool for automatically cleaning up gene alignments by \n\
        trimming to ORG and identifying and/or removing problematic sequences",
    help_template = "{before-help}\
        \n\nTool: {name}\
        \nDescription: {about-with-newline}\
        Version: {version}\
        \nAuthors: {author-with-newline}\
        \n{usage-heading} {usage}\
        \n{all-args}\
        {after-help}"
)]
pub struct Config {
    #[arg(short, long, help = "Input Fasta file",
    value_parser = validate_filename)]
    inp_fasta: String,

    #[arg(short, long, help = "Output Fasta file", 
        default_value_t = String::from("./output.fasta"),
    value_parser = validate_out_fasta)]
    out_fasta: String,

    #[arg(short, long, help = "The selected module(s)",
    value_parser = validate_modules)]
    module: String,
}

impl Config {
    pub(crate) fn get_module(&self) -> &str {
        &self.module
    }

    pub(crate) fn get_inp_fasta(&self) -> &str {
        &self.inp_fasta
    }

    pub(crate) fn get_out_fasta(&self) -> &str {
        &self.out_fasta
    }
}

fn validate_filename(name: &str) -> Result<String, String> {
    if name.is_empty() {
        Err(String::from("Filename cannot be empty"))
    } else if !Path::new(name).exists() {
        Err(String::from("Filename does not exist"))
    } else {
        Ok(name.to_string())
    }
}

fn validate_modules(module: &str) -> Result<String, String> {
    if module.is_empty() {
        Err(String::from("Module name cannot be empty"))
    } else if !["TrimToORF", "placeholder"].contains(&module) {
        Err(format!("'{}' not a known module", module))
    } else {
        Ok(module.to_string())
    }
}

fn validate_out_fasta(name: &str) -> Result<String, String> {
    if name.is_empty() {
        Err(String::from("Filename cannot be empty"))
    } else {
        Ok(name.to_string())
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn good_filename() {
        let test_name = "test_data/a_ha_h3_raw_500.fna";
        let result = validate_filename(test_name);
        assert_eq!(result, Ok("test_data/a_ha_h3_raw_500.fna".to_string()));
    }

    #[test]
    fn bad_filename() {
        let test_name = "notreal.fna";
        let result = validate_filename(test_name);
        assert_eq!(result, Err(String::from("Filename does not exist")));
    }

    #[test]
    fn no_filename() {
        let test_name = "";
        let result = validate_filename(test_name);
        assert_eq!(result, Err(String::from("Filename cannot be empty")));
    }

    #[test]
    fn good_module() {
        let test_module = "TrimToORF";
        let result = validate_modules(test_module);
        assert_eq!(result, Ok(test_module.to_string()));
    }

    #[test]
    fn bad_module() {
        let test_module = "NotAModule";
        let result = validate_modules(test_module);
        assert_eq!(
            result,
            Err(format!("'{}' not a known module", test_module))
        );
    }

    #[test]
    fn no_module() {
        let test_module = "";
        let result = validate_modules(test_module);
        assert_eq!(result, Err(String::from("Module name cannot be empty")));
    }
}
