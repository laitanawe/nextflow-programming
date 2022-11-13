---
title: "Nextflow DSL2 Subworkflow Definition"
teaching: 50
exercises: 10
questions:
- "What are the main features of a Nextflow module?"
- "What are the main components of a Nextflow pipeline?"
- "How do I run a Nextflow pipeline?"
objectives:
- "Explain the benefits of using Nextflow to develop pipelines."
- "Explain the components of a Nextflow process, pipeline."
- "Run a Nextflow pipeline."
keypoints:
- "A workflow is a sequence of tasks that process a set of data."
- "Nextflow scripts comprise of channels for controlling inputs and outputs, and processes for defining workflow tasks."
- "You run Nextflow modules and pipelines using the `nextflow run` command."
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


This is a Nextflow pipeline performing the following steps:
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

include { index } from './modules/index' // './modules/some_module' e.g. Index
include { quant } from './modules/quant' // './modules/other_module' e.g. Quant
include { fastqc } from './modules/fastqc' // './modules/another_more_module' e.g. FastQC

//  The default workflow
workflow rnaseq_quant {

take:
  transcriptome
  read_pairs_ch

main:
  index( transcriptome )
  fastqc( read_pairs_ch )
  quant( index.out, read_pairs_ch )

emit:
  quant.out
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
> > Launching `main.nf` [fervent_babbage] - revision: c54a707593
> > executor >  local (1)
> > [21/b259be] process > NUM_LINES (1) [100%] 1 of 1 ✔
> >
> >  ref1_1.fq.gz 58708
> > ~~~
> > {: .output}
> >
> > 1. The first line shows the Nextflow version number.
> > 1. The second line shows the run name `fervent_babbage` (adjective and scientist name) and revision id `c54a707593`.
> > 1. The third line tells you the process has been executed locally (`executor >  local`).
> > 1. The next line shows the process id `21/b259be`, process name, number of cpus, percentage task completion, and how many instances of the process have been run.
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
