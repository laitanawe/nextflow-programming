---
title: "Nextflow DSL2 Subworkflow Definition"
teaching: 50
exercises: 10
questions:
- "What are the main features of a subworkflowin Nextflow DSL2?"
- "What are the main components of a subworkflow in Nextflow DSL2?"
- "How do I run a subworkflow in Nextflow DSL2?"
objectives:
- "Explain the benefits of using subworkflows to develop pipelines in Nextflow DSL2."
- "Explain the components of a subworkflow, pipeline in Nextflow DSL2."
- "Run a subworkflow in Nextflow DSL2."
keypoints:
- "A subworkflow is a sequence of tasks that process a set of data."
- "Nextflow DSL2 scripts comprise of channels for controlling inputs and outputs, and processes for defining workflow tasks."
- "You run Nextflow subworkflows and modules and pipelines using the `nextflow run` command."
---

## Our example script

We are now going to look at a sample Nextflow script that performs some RNA-seq tasks.

For the Nextflow course demo, you need to download some files to follow the lesson.
1. Download <a href="https://laitanawe.github.io/nextflow-novice/data/nextflow-nov-lesson.tar.gz">nextflow-nov-lesson.tar.gz</a>
2. Make a directory called "nfdemo" on the Desktop
3. Move the nextflow-nov-lesson.tar.gz inside the nfdemo directory and cd into the nfdemo directory.
4. Unzip/extract `nextflow-nov-lesson.tar.gz` by typing `tar -xzvf nextflow-nov-lesson.tar.gz`
You should end up certain files within the folder **`nfdemo/data/ggal`** on your Desktop.


Open the file `main.nf` in the script directory with your favourite text editor.

This is a Nextflow script. It contains;

1. An optional interpreter directive ("Shebang") line, specifying the location of the Nextflow interpreter.
1. `nextflow.enable.dsl=2` to enable DSL2 syntax.
1. A multi-line Nextflow comment, written using C style block comments,
followed by a single line comment.
1. A pipeline parameter `params.reads` which is given a default value, of the relative path to the location of paired-end fastq files, as a string.
1. An unnamed `workflow` execution block, which is the default workflow to run.
1. A Nextflow channel used to read in data to the workflow.
1. A call to the process `index`.
1. An operation on the process output, using the channel operator `view()`.
1. Nextflow `process` blocks named `fastqc`, `quant`, `multiqc`, which defines what the process does.
1. An `input` definition block that assigns the input to a variable, and declares that it should be interpreted as a file `path`.
1. An `output` definition block that uses the path object (we can also use Linux/Unix standard output stream `stdout`) from the script block.
1. A `script` block that contains some bash commands.


This is a Nextflow subworkflow (./rnaseq_mods.nf) performing the following steps:
1. Indexes a transcriptome file.
1. Performs quality control.
1. Performs quantification.
1. Creates a MultiQC report.

~~~
// Save your ./rnaseq_mods.nf script so that it looks like this:

    include { index } from './modules/index' // './modules/some_module' //Cut and Paste index block into ./modules/index.nf
    include { quant } from './modules/quant' // './modules/other_module' //Cut and Paste index block into ./modules/quant.nf
    include { fastqc } from './modules/fastqc' // './modules/another_more_module' //Cut and Paste index block into ./modules/fastqc.nf
    include { multiqc } from './modules/multiqc' // './modules/one_more_module' //Cut and Paste index block into ./modules/multiqc.nf

//  The default subworkflow
workflow rnaseq_sub {

take:
  transcriptome
  read_pairs_ch

main:
  index( transcriptome )
  fastqc( read_pairs_ch )
  quant( index.out, read_pairs_ch )

emit:
  fastqc.out.mix( quant.out ).collect()
}


// Copy from the shebang line to the end so that your main.nf script now looks like this:

#!/usr/bin/env nextflow

/* Contributors:
 * - Awe, O.I
 */

nextflow.enable.dsl = 2

/*
 * Default pipeline parameters can be overriden on the command line eg.
 * given 'params.foo' specify on the run command line '--foo some_value'
 */

/* This is a general pattern for pairs of files. Triplets are also possible.
 * To point to other samples apart from ggal_gut, you can use a glob pattern: "$baseDir/data/ggal/*_{1,2}.fq"
 */
params.reads = "$baseDir/data/ggal/ggal_gut_{1,2}.fq"
params.transcriptome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.outdir = "results"
params.multiqc = "$baseDir/multiqc"

/* log.info is the same thing as println but Nextflow also puts it in the log */

log.info """\
RNASEQ NEXTFLOW PIPELINE uses the ffg bioinformatics tools: Salmon, FastQC, MultiQC
transcriptome : ${params.transcriptome}
reads         : ${params.reads}
outdir        : ${params.outdir}
"""

/*  Comments are uninterpreted text included with the script.
    They are useful for describing complex parts of the workflow
    or providing useful information such as workflow usage.

    Usage:
       nextflow run wc.nf --input <input_file>

    Multi-line comments start with a slash asterisk /* and finish with an asterisk slash. */
//  Single line comments start with a double slash // and finish on the same line

/*  Workflow parameters are written as params.<parameter>
    and can be initialised using the `=` operator. */
// There is a docker container with salmon in: https://hub.docker.com/r/nextflow/rnaseq-nf
// To execute by using docker by default and skip -with-docker at the cmd line, add this to your config: docker.enabled = true
// In your nextflow.config file: process.container = 'nextflow/rnaseq-nf'
// Or at the top of your process definition/script, container = 'nextflow/rnaseq-nf'
// docker.runOptions='-u $(id -u):$(id -g)'
// To run your salmon using the docker container, you can do: nextflow run script.nf -with-docker
// To use a docker image different from the one specified in the config, you can do: nextflow run script.nf -with-docker repo_name/image:tag
// For execution caching / resume, you can do: nextflow run script.nf -with-docker repo_name/image:tag -resume
// For more cpus, you can do: nextflow run script.nf -with-docker repo_name/image:tag -process.cpus 4

/*  A Nextflow include line for the subworkflow block
    Process names are written, by convention, in uppercase but it could also be in lowercase.
    This convention is used to enhance workflow readability. */

    include { rnaseq_sub; multiqc } from './rnaseq_mods' // './a_subworkflow' //Cut and Paste index,fastqc,quant,multiqc blocks into ./rnaseq_mods.nf

//  The default workflow
workflow {
// Input channel objects can be specified here in the workflow scope.
// Input data is received through channels
// input_ch = Channel.fromPath(params.input)
read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists: true )
// read_pairs_ch channel is a tuple where the 1st item is a value or the sampleID and the 2nd tuple item is a list of paths
read_pairs_ch.view()

// This section is an example on how to convert a tuple into a list, if you're interested.
Channel.fromFilePairs( params.reads, flat: true )
// read_pairs_flat_ch channel can be turned into a flat list
       .set{read_pairs_flat_ch}
read_pairs_flat_ch.view()

/*  The script to execute is called by its process name,
    and input is provided between brackets. */

/*  Process output is accessed using the `out` channel.
    The channel operator view() can be used to print a
    process' output to the terminal.
    The operator, .mix() combines all items from multiple channels and the resulting number of items is an addition of the channels.
    The operator, .collect() flattens a list and converts it into a single element/item leading to one task in the process.
     */

    rnaseq_sub( params.transcriptome, read_pairs_ch )
    multiqc( rnaseq_sub.out )

}

~~~~
{: .language-groovy}

To run a Nextflow script use the command `nextflow run <script_name>`.

> ## Run a Nextflow  script
> Run the script by entering the following command in your terminal:
>
> ~~~
> $ nextflow run main.nf
> ~~~
> {: .language-bash}
> > ## Solution
> > You should see output similar to the text shown below:
> >
> > ~~~
> > N E X T F L O W  ~  version 20.10.0
> > Launching `main.nf` [sleepy_avogadro] DSL2 - revision: 83e5d597be
> > RNASEQ NEXTFLOW PIPELINE uses the ffg bioinformatics tools: Salmon, FastQC, MultiQC
> > transcriptome : /home/aweo/Desktop/nfdemo/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa
> > reads         : /home/aweo/Desktop/nfdemo/data/ggal/ggal_gut_{1,2}.fq
> > outdir        : results
> > executor >  local (4)
> > [51/a6b821] process > rnaseq_sub:index (ggal_1_48850000_49020000) [100%] 1 of 1 ✔
> > [0c/1f6c84] process > rnaseq_sub:fastqc (FASTQC on ggal_gut)      [100%] 1 of 1 ✔
> > [46/a7019e] process > rnaseq_sub:quant (1)                        [100%] 1 of 1 ✔
> > [7c/2f055d] process > multiqc                                     [100%] 1 of 1 ✔
> > [ggal_gut, [/home/aweo/Desktop/nfdemo/data/ggal/ggal_gut_1.fq, /home/aweo/Desktop/nfdemo/data/ggal/ggal_gut_2.fq]]
> > [ggal_gut, /home/aweo/Desktop/nfdemo/data/ggal/ggal_gut_1.fq, /home/aweo/Desktop/nfdemo/data/ggal/ggal_gut_2.fq]
> > ~~~
> > {: .output}
> >
> > 1. The first line shows the Nextflow version number.
> > 1. The second line shows the run name `fervent_babbage` (adjective and scientist name) and revision id `c54a707593`.
> > 1. The third line tells you the process has been executed locally (`executor >  local`).
> > 1. The next line shows the log info followed by the process ids e.g. `51/a6b821`, process name, number of cpus, percentage task completion, and how many instances of the process have been run. Take note of subworkflow process names rnaseq_sub:index, rnaseq_sub:fastqc, and rnaseq_sub:quant.
> > 1. The final line is the output of the `view` operator.
> {: .solution}
{: .challenge}


> ## Process identification
> The hexadecimal numbers, like 61/1f3ef4, identify the unique process execution.
> These numbers are also the prefix of the directories where each process is executed.
> You can inspect the files produced by changing to the directory `$PWD/work` and
> using these numbers to find the process-specific execution path. We will learn how exactly
> nextflow using *work* directory to execute processes in the following sections.
{: .callout}



{% include links.md %}
