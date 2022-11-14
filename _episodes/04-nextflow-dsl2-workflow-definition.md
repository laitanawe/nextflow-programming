---
title: "Nextflow DSL Workflow Definition"
teaching: 50
exercises: 10
questions:
- "What are the main differences between a Nextflow DSL2 and DSL1 process?"
- "What are the main components of a Nextflow DSL2 pipeline?"
- "How do I run a Nextflow DSL2 pipeline?"
objectives:
- "Explain the benefits of using Nextflow DSL2 to develop pipelines."
- "Explain the components of a Nextflow DSL2 process, pipeline."
- "Run a Nextflow DSL2 pipeline."
keypoints:
- "A DSL2 workflow is a sequence of tasks that processes a set of data."
- "In Nextflow DSL2 processes, the process definition is no longer tied to specific channels (as in DSL1) and so, a DSL2 module can be used independently. This enables reuse of DSL2 components. From some module, we can include DSL2 processes in a workflow."
- "You run Nextflow DSL2 modules and pipelines using the `nextflow run` command."
---

## Our example script

We are now going to look at a sample Nextflow DSL2 workflow that performs some RNA-seq tasks.

For the Nextflow course demo, you need to download some files to follow the lesson.
1. Download <a href="https://laitanawe.github.io/nextflow-novice/data/nextflow-nov-lesson.tar.gz">nextflow-nov-lesson.tar.gz</a>
2. Make a directory called "nfdemo" on the Desktop
3. Move the nextflow-nov-lesson.tar.gz inside the nfdemo directory and cd into the nfdemo directory.
4. Unzip/extract `nextflow-nov-lesson.tar.gz` by typing `tar -xzvf nextflow-nov-lesson.tar.gz`
You should end up certain files within the folder **`nfdemo/data/ggal`** on your Desktop.

Open the file `main.nf` in the script directory with your favourite text editor.

This is a Nextflow DSL2 pipeline performing the following steps:
1. Indexes a transcriptome file.
1. Performs quality control.
1. Performs quantification.
1. Creates a MultiQC report.

~~~
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

    Multi-line comments start with a slash asterisk /* and finish with an asterisk slash.
    */
//  Single line comments start with a double slash // and finish on the same line

/*  Workflow parameters are written as params.<parameter>
    and can be initialised using the `=` operator.
    */
// There is a docker container with salmon in: https://hub.docker.com/r/nextflow/rnaseq-nf
// To execute by using docker by default and skip -with-docker at the cmd line, add this to your config: docker.enabled = true
// In your nextflow.config file: process.container = 'nextflow/rnaseq-nf'
// Or at the top of your process definition/script, container = 'nextflow/rnaseq-nf'
// docker.runOptions='-u $(id -u):$(id -g)'
// To run your salmon using the docker container, you can do: nextflow run script.nf -with-docker
// To use a docker image different from the one specified in the config, you can do: nextflow run script.nf -with-docker repo_name/image:tag
// For execution caching / resume, you can do: nextflow run script.nf -with-docker repo_name/image:tag -resume
// For more cpus, you can do: nextflow run script.nf -with-docker repo_name/image:tag -process.cpus 4

/* Module: a single atomic process or Nf script containing 1+ processes that can be imported from another Nf script.
   Sub-workflow: a few chained modules or Nf script containing 1+ processes that can be imported from another Nf script.
   Workflow: an end-to-end pipeline.
   A process/sub-workflow/workflow can be a module.
*/
/*  A Nextflow process block, using include statements:
    Curly brackets are mandatory in include statements, so you can include specific/multiple processes from a file e.g.
    include { index; quant } from './modules/quant'
    from ./some_module or from ./some_module.nf or from /some/module
    The .nf is optional but the path to the module file has to be relative or absolute
    */

    include { index } from './modules/index' // './modules/some_module' //Cut and Paste index block into ./modules/index.nf
    include { quant } from './modules/quant' // './modules/other_module' //Cut and Paste index block into ./modules/quant.nf
    include { fastqc } from './modules/fastqc' // './modules/another_more_module' //Cut and Paste index block into ./modules/fastqc.nf
    include { multiqc } from './modules/multiqc' // './modules/one_more_module' //Cut and Paste index block into ./modules/multiqc.nf

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

    index( params.transcriptome )
    fastqc( read_pairs_ch )
    quant( index.out, read_pairs_ch )
    multiqc( fastqc.out.mix( quant.out ).collect() )

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
> > [3b/a13c45] process > index (ggal_1_48850000_49020000) [100%] 1 of 1 ✔
> > [ff/f60afe] process > fastqc (FASTQC on ggal_gut)      [100%] 1 of 1 ✔
> > [34/84f09d] process > quant (1)                        [100%] 1 of 1 ✔
> > [aa/25f8e4] process > multiqc                          [100%] 1 of 1 ✔
> > [ggal_gut, [/home/aweo/Desktop/nfdemo/data/ggal/ggal_gut_1.fq, /home/aweo/Desktop/nfdemo/data/ggal/ggal_gut_2.fq]]
> > [ggal_gut, /home/aweo/Desktop/nfdemo/data/ggal/ggal_gut_1.fq, /home/aweo/Desktop/nfdemo/data/ggal/ggal_gut_2.fq]
> > ~~~
> > {: .output}
> >
> > 1. The first line shows the Nextflow version number.
> > 1. The second line shows the run name `fervent_babbage` (adjective and scientist name) and revision id `c54a707593`.
> > 1. The third line tells you the process has been executed locally (`executor >  local`).
> > 1. The next line shows the log info followed by the process ids e.g. `3b/a13c45`, process name, number of cpus, percentage task completion, and how many instances of the process have been run.
> > 1. The final line is the output of the `view` operator.
> {: .solution}
{: .challenge}


> ## Process identification
> The hexadecimal numbers, like 3b/a13c45, identify the unique process execution.
> These numbers are also the prefix of the directories where each process is executed.
> You can inspect the files produced by changing to the directory `$PWD/work` and
> using these numbers to find the process-specific execution path. We will learn how exactly
> nextflow using *work* directory to execute processes in the following sections.
{: .callout}



{% include links.md %}
