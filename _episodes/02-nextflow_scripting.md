---
title: "Nextflow scripting"
teaching: 40
exercises: 20
questions:
- "What language are Nextflow scripts written in?"
- "How do I store values in a Nextflow script?"
- "How do I write comments in a Nextflow script?"
- "How can I store and retrieve multiple values?"
- "How are strings evaluated in Nextflow?"
- "How can I create simple re-useable code blocks?"
objectives:
- "Understand what language Nextflow scripts are written in."
- "Define variables in a script."
- "Create lists of simple values."
- "Comment Nextflow scripts."
- "Explain what a list is."
- "Explain what string interpolation is."
- "Understand what a closure is."

keypoints:
- "Nextflow is a Domain Specific Language (DSL) implemented on top of the Groovy programming language."
- "To define a variable, assign a value to it e.g., `a = 1` ."
- "Comments use the same syntax as in the C-family programming languages: `//` or multiline `/* */`. "
- "Multiple values can be stored in lists [value1, value2, value3, ...] or maps [chromosome: 1, start :1]."
- "Lists are indexed and sliced with square brackets (e.g., list[0] and list[2..9])"
- "String interpolation (variable interpolation, variable substitution, or variable expansion) is the process of evaluating a string literal containing one or more placeholders, yielding a result in which the placeholders are replaced with their corresponding values."
- "A closure is an expression (block of code) encased in `{}` e.g. `{ it * it }`."

---

Nextflow is a Domain Specific Language (DSL) implemented on top of the Groovy programming language, which in turn is a super-set of the Java programming language. This means that Nextflow can run any Groovy and Java code. It is not necessary to learn Groovy to use Nextflow DSL but it can be useful in edge cases where you need more functionality than the DSL provides.

## Language Basics

### Printing values

To print something is as easy as using the `println` method (`println` is a compression of "print line") and passing the text to print in quotes.
The text is referred to as a `string` as in a string of characters.
~~~
println("Our script works!")
~~~
{: .language-groovy }

~~~
Our script works!
~~~
{: .output }

**Parenthesis** for function invocations are optional. Therefore also the following is a valid syntax.

~~~
println "Our script works without parenthesis!"
~~~
{: .language-groovy }

~~~
Our script works without parenthesis!
~~~
{: .output }

## Methods

`println` is an example of a Groovy method. A method is just a block of code which only runs when it is called.
You can pass data, known as parameters, into a method using the method name followed by brackets `()`.
Methods are used to perform certain actions, and they are also known as functions.
Methods enable us to reuse code: define the code once, and use it many times.

## Comments

When we write any code it is useful to document it using comments.
In Nextflow comments use the same syntax as in the C-family programming languages.
This can be confusing for people familiar with the `#` syntax for commenting in other languages.

~~~
// This is a single line comment. Everything after the // is ignored.

/*
   Comments can also
   span multiple
   lines.
 */
~~~
{: .language-groovy }


## Variables

In any programming language, you need to use variables to store different types of information. A variable is a pointer to a space in the computer's memory that stores the value associated with it.

Variables are assigned using `=` and can have any value. Groovy is dynamically-typed which means the variable's data type is based on its value. For example, setting `x = 1` means `x` is an integer number, but if it is later set to `x = "hello"` then it becomes a String.

> ## Variable scope
> When we create a variable using the `x = 1` syntax we can access, (`scope`), it anywhere (`globally`) in the script. A variable declared in this fashion is sometimes called a public variable.
>
> We can also define variables with a data `type` e.g. `String x="Hello"` or with the `def` keyword `def x=1`. This effects the accessibility (`scope`) of the variable.
> This is called lexical scoping (sometimes known as static scoping) that sets the scope  of a variable so that it may only be accessed from within the block of code in which it is defined. A variable declared in this fashion is sometimes called a private variable.
{: .callout }

### Types of Data

Groovy knows various types of data. four common ones are:


* `String` − These are text literals which are represented in the form of chain of characters enclosed in quotes. For example `"Hello World"`.
* `int` − This is used to represent whole numbers. An example is `1234`.
* `Boolean` − This represents a Boolean value which can either be `true` or `false`.
* `float` - This is used to represent floating point number `12.34` .

A more complete list can be found [here](https://www.tutorialspoint.com/groovy/groovy_data_types.htm)


In the example below, variable `my_var` has an integer value of `1`:

~~~
//int − This is used to represent whole numbers.
my_var = 1
~~~
{: .language-groovy }

To create a variable with a floating point value, we can execute:

~~~
//float − This is used to represent floating point numbers.
my_var = 3.1499392
~~~
{: .language-groovy }

To create a Boolean value we assign the value `true` or `false`.  
**Note:* Do not enclose a Boolean value in quotes or they will be interpreted as a string.

~~~
//Boolean − This represents a Boolean value which can either be true or false.
my_var = false
~~~
{: .language-groovy }


And to create a string, we add single or double quotes around some text.

For example:

~~~
//String - These are text literals which are represented in the form of chain of characters
my_var = "chr1"
~~~
{: .language-groovy }

### Multi-line strings

A block of text that span multiple lines can be defined by delimiting it with triple single `'''` or double quotes `"""`:

~~~
text = """
    This is a multi-line string
    using triple quotes.
    """
~~~
{: .language-groovy }


To display the value of a variable to the screen in Groovy, we can use the `println` method passing the variable name are a parameter.

~~~
x = 1
println(x)
~~~
{: .language-groovy }


~~~
1
~~~
{: .output }

### Slashy strings

Strings can also be defined using the forward slash `/` character as delimiter. They are known as `slashy strings` and are useful for defining regular expressions and patterns, as there is no need to escape backslashes e.g `\n` specifies a new line. As with double quote strings they allow to interpolate variables prefixed with a `$` character.

Try the following to see the difference:

~~~
x = /ATP1B2\TP53\WRAP53/
println(x)
~~~
{: .language-groovy }

~~~
ATP1B2\TP53\WRAP53
~~~
{: .output }

~~~
y = 'ATP1B2\TP53\WRAP53'
println(y)
~~~
{: .language-groovy }

Produces an error as the `\` is a special characters that we need to escape.

~~~
// use \ to escape
y = 'ATP1B2\\TP53\\WRAP53'
println(y)
~~~
{: .language-groovy }

~~~
ATP1B2\TP53\WRAP53
~~~
{: .output }

### String interpolation

To use a variable inside a single or multi-line double quoted string `""`  prefix the variable name with a `$` to show it should be interpolated:

~~~
chr = "1"

//When using your variable within double quotes, you must prefix the $ for your variable to interpolate
println("processing chromosome $chr")

//You don't need the $ prefix before chr because it's unquoted when concatinating as seen below:              
println("processing chromosome " + chr)

The examples below highlight how to interpolate in Nextflow:

var_int = 17

//No $ prefix but there is interpolation:
println(var_int)

//Double Quotes, No $, No Interpolation:
println("var_int")

//Double Quotes, $ Interpolation:
println("$var_int")

//Double Quotes, ${} Interpolation:
println("${var_int}")

//Double Quotes, No $, No Interpolation:
println("{var_int}")

//Single Quotes, ${} but No Interpolation:
println('${var_int}')

~~~
{: .language-groovy }

~~~
processing chromosome 1

17
var_int
17
17
{var_int}
${var_int}

~~~
{: .output }

**Note:** Variable names inside single quoted strings do not support String interpolation.

~~~
chr = "1"
//Your variable will not interpolate when used within single quotes
println('processing chromosome $chr')

~~~
{: .language-groovy }

~~~
processing chromosome $chr
~~~
{: .output }

### String interpolation within the script block

In the script block, to use a nextflow variable inside a single or multi-line double quoted string `""` prefix the variable name with a `$` to show it should be interpolated.

<b>NOTE: a.) Within script block</b>:
`$` is <b>mandatory for interpolation within the script block and `""` or `''` is optional when you prefix a Nextflow variable with `$` in the script block</b>. <b>Even when using `$` within a bash `#` comment</b> in the script block, you must use `\$` in order to treat the `$` as part of the comment.

<b>NOTE: b.) Outside the script block</b>:
Outside the script block and in the Nextflow scope (<i>any area that is not within the script block</i>), `""` is required if you want to interpolate a Nextflow variable using a `$` prefix. Also, in the Nextflow scope, use a Nextflow variable name within the arguments of a groovy function/operator if there's no `""` i.e. you don't need a `$` prefix for interpolation in that case.

~~~
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process testme {

input:
 val mystr

output:
 stdout

script:
 """
 #The examples below highlight how to interpolate in the script block:

 echo mystr #No \$ prefix, and there is No interpolation.

 echo "mystr" #Double Quotes, No \$, No Interpolation.

 echo "$mystr" #Double Quotes, \$ Interpolation.

 echo "${mystr}" #Double Quotes, \${} Interpolation.

 echo "{mystr}" #Double Quotes, No \$, No Interpolation.

 echo '${mystr}' #Single Quotes, \${} but there is Interpolation.

 echo '$mystr' #Single Quotes, \$ but there is Interpolation.

 echo ${mystr} #No Quotes, \${} but there is Interpolation.

 """
}

workflow {

any_var_name = "Our Script Works!"
params.input = Channel.from(any_var_name)
testme(params.input)
testme.out.view()

}

~~~
{: .language-groovy }

~~~

executor >  local (1)
[86/2f14f5] process > testme (1) [100%] 1 of 1 ✔
mystr
mystr
Our Script Works!
Our Script Works!
{mystr}
Our Script Works!
Our Script Works!
Our Script Works!

~~~
{: .output }

**Note:** Variable names inside single quoted strings do not support String interpolation.

## Lists

To store multiple values in a variable we can use a List.
A List  (also known as array) object can be defined by placing the list items in square brackets and separating items by commas `,`:

~~~
kmers = [11,21,27,31]
~~~
{: .language-groovy }

You can access a given item in the list with square-bracket notation `[]`. These positions are numbered starting at 0, so the first element has an index of 0.

~~~
kmers = [11,21,27,31]
println(kmers[0])

~~~
{: .language-groovy }

~~~
11
~~~
{: .output}

We can use negative numbers as indices in Groovy. They count from the end of the list rather than the front: the index `-1` gives us the last element in the list, `-2` the second to last, and so on. Because of this, `kmers[3]` and `kmers[-1]` point to the same element in our example list.

~~~
#!/usr/bin/env nextflow
kmers = [11,21,27,31]
//Lists can also be indexed with negative indexes
println(kmers[3])
println(kmers[-1])

println("The last element of kmers is: " + kmers[3])

println("The last element of kmers is: $kmers[-1]")

println("The last element of kmers is: " + kmers[0..2])

println("The last element of kmers is: ${kmers[-1]}")

println("The second-to-the-last element of kmers is: ${kmers[-2]}")

~~~
{: .language-groovy }
~~~
N E X T F L O W  ~  version 22.10.1
Launching `script5.nf` [grave_hodgkin] DSL2 - revision: 249ba9a9f4
31
31
The last element of kmers is: 31
The last element of kmers is: [11, 21, 27, 31][-1]
The last element of kmers is: [11, 21, 27]
The last element of kmers is: 31
The second-to-the-last element of kmers is: 27

~~~
{: .output}

Lists can also be indexed using a range. A range is a quick way of declaring a list of consecutive sequential numbers.
To define a range use `<num1>..<num2>` notation.

~~~
kmers = [11,21,27,31]
// The first three elements using a range.
println(kmers[0..2])
// Because kmers is a list, this will print the first three elements of the kmers list
~~~
{: .language-groovy }
~~~
[11, 21, 27]
~~~
{: .output}

### String interpolation of list elements

To use an expression like `kmer[0..2]` inside a double-quoted string `""`, we use the `${expression}` syntax, similar to Bash shell scripts.

For example, the expression below without the `{}`""

~~~
kmers = [11,21,27,31]
println("The first three elements in the Lists are: $kmers[0..2]")
// Because $kmers is a list, this will print the kmers list, followed by [0..2]
~~~
{: .language-groovy }

would output.

~~~
The first three elements in the Lists are: [11, 21, 27, 31][0..2]
~~~
{: .output}

We need to enclose the `kmers[0..2]` expression inside `${}` as below to get the correct output.

~~~
kmers = [11,21,27,31]
println("The first three elements in the Lists are: ${kmers[0..2]}")
~~~
{: .language-groovy }


~~~
The first three elements in the Lists are: [11, 21, 27]
~~~
{: .output}


### List methods

Lists have a number of useful methods that can perform operations on their contents. See more [here](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/List.html). When using a method on a type of object you need prefix the method with the variable name.

For example, in order to get the length of the list use the list `size` method:

~~~
mylist = [0,1,2]

println(mylist.size())

//inside a string need we need to use the ${} syntax
println("list size is: ${mylist.size()}")

// When there's no double quotes, you can use this:
println("list size is: " + mylist.size())

~~~
{: .language-groovy }

~~~
3
list size is: 3
~~~
{: .output }

We can use the `get` method items to retrieve items in a list.

~~~
mylist = [100,150,200]
println mylist.get(1)

println "${mylist.get(1)}"

println(mylist.get(1))

println("${mylist.get(1)}")

~~~
{: .language-groovy }

~~~
150
150
150
~~~
{: .output }

Listed below are a few more common list methods and their output on a simple example.

~~~
mylist = [1,2,3]
println mylist
println mylist + [1]
println mylist - [1]
println mylist * 2
println mylist.reverse()
println mylist.collect{ it+3 }
println mylist.unique().size()
println mylist.count(1)
println mylist.min()
println mylist.max()
println mylist.sum()
println mylist.sort()
println mylist.find{it%2 == 0}
println mylist.findAll{it%2 == 0}
~~~
{: .language-groovy }

~~~
[1, 2, 3]
[1, 2, 3, 1]
[2, 3]
[1, 2, 3, 1, 2, 3]
[3, 2, 1]
[4, 5, 6]
3
1
1
3
6
[1, 2, 3]
2
[2]
~~~
{: .output }

> ## Create List and retrieve value
> Create a list object `list` with the values 1 to 10.
> Access the fifth element in the list using with square-bracket notation or using the `get` method and
> print the results
> > ## Solution
> > ~~~
> > list = [1,2,3,4,5,6,7,8,9,10]
> > //or
> > list = 1..10
> > println("${list[4]}")
> > //or
> > println("${list.get(4)}")
> > ~~~
> > {: .language-groovy }
> > The fifth element is `5`. Remember that the array index starts at 0.
> > {: .output}
> {: .solution}
{: .challenge}


## More resources

The complete Groovy language documentation is available at this [link](http://groovy-lang.org/documentation.html#languagespecification).

A great resource to master Apache Groovy syntax is Groovy in [Action](https://www.manning.com/books/groovy-in-action-second-edition).
