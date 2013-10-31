Documentation for LingPy
========================

If you contribute to LingPy, you should document your code.
The full documentation for each script in LingPy is collected in the folder ```lingpy/doc/source/docu```. In this folder LingPy's internal structure should be reflected and all pieces of code should be rendered here in the documentation. 
Apart from providing documentation directly through the functions and modules you write, this is the most important step for consistent documentation. Below, I try to summarize, how documentation should be provided in more detail

## How to provide documentation

### First step: Documentation within the code

The first step for documentation is the documentation within the code. For examples, please have a look at this example from the Wordlist classe: https://github.com/lingpy/lingpy/blob/master/lingpy/basic/wordlist.py#L35.

### Second step: Documentation within in the Docu-Sources

Apart from the first step, the documentation needs to be repeated in the lingpy/doc/source/docu-Folder. Here, you simply select the folder where your script is located, create a script if it is not already there, and type in the relevant sphinx-commands that are needed to document the code properly. These commands are best to be learnt from the examples of existing commands, as, for example, the examples on Wordlists, which you find from this link: http://lingpy.org/docu/basic/generated/lingpy.basic.wordlist.Wordlist.html#lingpy.basic.wordlist.Wordlist. If you go to the right page, you can see the source by simply pressing "show source".

## Why provide documentation?

Documentation needs to be provided in order to enable other users to understand and use our code. Only those code-pieces which are substantially documented will be displayed on the documentation side on http://lingpy.org. So there's some good reason, if you want your code to be used, to put some effort into documenting it properly. 


