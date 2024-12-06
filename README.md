As a first step, this document is under governance review. When the review completes as appropriate per local and agency processes, the project team will be allowed to remove this notice. This material is draft.

# clean-genes overview

## Project goal
To automatically clean up a gene alignment by trimming to open reading frame (ORF) and identifying and/or removing problematic sequences. Basically an alignment normalization tool that expects genes as inputs and could have specific value for Influenza genes. This is intended to be a Rust crate designed to do a lot of work as quickly as possible.

## Project title
current proposal: clean-genes

Justification:
1. The name is available as a software in general and as a Rust crate on crates.io
2. Cleaning is what the software does
3. Genes are the expected input

## Motivation

I have done this work more/less manually before and it takes a lot of time, is hard to explain, and is very susceptible to human error. Where I did this work a lot was when I started working with influenza at USDA-ARS. Because influenza genomes are small and a deep study can involve a great number of viruses this work was frequently necessary

## Project plan

I intend to build the software in a modular fashion such that it can:
1. Be built in smaller pieces to demonstrate value iteratively
2. Be utilized by people with different needs from what I'm familiar with
3. Run more efficiently when the user does not want or need all available functionality

### Potential modules
- ORF trimming: Determines the reading frame based on stop-codon frequency and conservation. Optionally trims the sequence using these determined reading frames
- Frameshift ID: Identifies sequences with frameshift mutations based on stop-codon frequency and conservation. Optionally removes these sequences and prompts the user to re-align and run this software again.
- Poor Quality Filter: Identifies sequences with an abundance of missing data, overly short length, or significant gaps. Optionally removes these sequences
- Taxonomic Outlier ID: Identifies sequences that don't belong in the alignment taxonomically. Could potentially use GC content, k-mers, and/or sequence similarity. I am resistant to using reference sequences for this, but it's not impossible. Optionally removes these sequences.
- Codon Optimized ID: Identifies codon-optimized sequences like those created for vaccines based on codon distribution in the alignment and/or an understanding of what a normal codon distribution is. Optionally removes these sequences.
- Gap Groups: Identifies groups of sequences with shared indel patterns. Also identifies outliers. Optionally outputs these sequences in seperate fasta files excluding the common gaps. 

When used together, the removal steps are no longer optional. For ID-only the modules must be run one-at-a-time.

Alignment would likely be done in MAFFT. MAFFT is fairly fast and produces high-quality results consistently using defualt parameters.

### Phase 1
Build the basic user-interface and the ORF trimming module.
- ORF trimming

### Phase 2
Add the Poor Quality Filter module.

workflow if all are included:
1. Alignment (if needed)
2. Poor Quality Filter
3. Re-alignment
4. ORF trimming

### Phase 3
Add the Gap Groups module

workflow if all are included:
1. Alignment (if needed)
2. Poor Quality Filter
3. Re-alignment
4. Gap Groups
5. ORF trimming

### Phase 4
Add the Frameshift ID module

1. Alignment (if needed)
2. Poor Quality Filter
3. Re-Alignment
4. Gap Groups
5. ORF trimming
6. Frameshift ID
7. Re-alignment

### Phase 5
Add the Codon Optimized ID module

workflow if all are included:
1. Alignment (if needed)
2. Poor Quality Filter
3. Re-Alignment
4. Gap Groups
5. ORF trimming
6. Frameshift ID
7. Re-alignment
8. Codon Optimized ID
9. Re-alignment

### Phase 6
Add the Taxonomic Outlier module

workflow if all are included:
1. Alignment (if needed)
2. Poor Quality Filter
3. Re-Alignment
4. Gap Groups
5. ORF trimming
6. Frameshift ID
7. Re-alignment
8. Codon Optimized ID
9. Re-alignment
10. Taxonomic Outlier 
11. Re-alignment


### Timeline
The time required for software development is very hard to predict. Things like VCM, holidays, unforseen challenges, and conferences would delay the schedule. Given that I am not solely focused on this project and it is holiday season I believe I can complete Phase 1 by the end of December. Optimization for performance may or may not be done before Phase 2, but correctness will be ensured. 

Phase 2 will likely be a little harder than Phase 1 due to multiple tasks being completed in 1 module as well as integration of alignment and could take a couple of months.

Phase 3 should be similar in complexity to Phase 1 taking about a month.

Phase 4 is seems like a bit of a harder problem and could take a couple of months.

Phases 5 and 6 will be very challenging and I am hesitant to put a timeline on them at all.

Unit tests, Regression tests, benchmarking, and fuzzing will all be used to ensure correctness and maximize performance.

Between each phase I intend to discuss my progress with my supervisor Brian Lee, my Rust mentor Sam Shepard, and other interested colleauges including Allen Kim who has previously expressed interest in contributing to this project. 

## Notices

### Public Domain Standard Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

### License Standard Notice
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

### Privacy Standard Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](DISCLAIMER.md)
and [Code of Conduct](code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

### Contributing Standard Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

### Records Management Standard Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

### Additional Standard Notices
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template) for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/main/CONTRIBUTING.md), [public domain notices and disclaimers](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md), and [code of conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md).
