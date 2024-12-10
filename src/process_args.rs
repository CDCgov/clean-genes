use clap::Parser;

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
    value_parser = validate_filename)]
    out_fasta: String,

    #[arg(short, long, help = "The selected module(s)",
    value_parser = validate_modules)]
    module: String,
}

impl Config {
    pub(crate) fn get_module(&self) -> &str {
        &self.module
    }
}

fn validate_filename(name: &str) -> Result<String, String> {
    if name.is_empty() {
        Err(String::from("Filename cannot be empty"))
    } else {
        Ok(name.to_string())
    }
}

fn validate_modules(module: &str) -> Result<String, String> {
    if module.is_empty() {
        Err(String::from("Module name cannot be empty"))
    } else if !["TrimToORF", "placeholder"].contains(&module) {
        Err(format!("'{}' not a known module", module))
    } else {
        Ok(module.to_string())
    }
}
