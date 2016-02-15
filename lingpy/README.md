# LingPy

## Basic Modules

The library in it's current state consist of the following modules:

* algorithm - directory dedicated to specific algorithms which are mainly used program-internally
* align - basic module for the conduction of alignment analyses
* basic - basic module in which general classes are defined that are crucial for most analyses
* compare - basic module for the comparison of languages, dialects, doculects
* data - basic module for the handling of externally defined data (specifications of IPA, sound-class models, configuration files, name-spaces etc.)
* evaluate - module for the evaluation of algorithms (alignment, automatic cognate judgments)
* meaning - module is dedicated to semantic approaches to automatic language comparison
* read - basic module defines specific functions to read in various data types
* sequence - basic module for the handling of sequential data (most often: words)

## Global Variables

When loading LingPy in a python session using either the 
```python
>>> from lingpy import *
```
or the 
```python
>>> import lingpy as lp
```
statements, not only functions are loaded, but also variables containing static data. For example, 
```python
>>> asjp
<sca-model "asjp">
```
refers to the sound-class model "asjp". This model consists of a so-called "converter", a dictionary that converts IPA to an internal sound-class alphabet consisting of quasi-asjp characters (close to the characters of the [ASJP project](http://email.eva.mpg.de/~wichmann/ASJPHomePage.htm])). The global variables defined in such a way are the following:

* asjp: ASJP sound-class model
* sca: SCA sound-class model
* dolgo: Dolgopolsky's original sound-class model (slightly adapted to LingPy's requirements)
* ipa_vowels: All vowel characters recognized by LingPy as characters referring to vowels in the IPA
* ipa_diacritics: All characters which LingPy recognizes as IPA diacritics
* ipa_tones: All standard tone characters, including sub- and superscript numerals.
* `_color`: A color scheme, that converts IPA characters to 10 specific colors as defined by the Dolgopolsky sound-class model

As an alternative, there is a dictionary called "rcParams", containing all these variables, and additional ones. This dictionary is used for internal coding purposes and stores parameters that are globally set (if not defined otherwise by the user), such as

* specific debugging messages (warnings, etc.)
* specific flags (verbose, debug)
* default values, such as "gop" (gap opening penalty for alignment analyses), scale (scaling factor by which extended gaps are penalized), or "figsize" (the default size of figures if data is plotted using matplotlib

These default values can be changed with help of the ```rc``` function that takes any keyword and any variable as input and adds or modifies the specific key of the rcParams dictionary, but also provides more complex functions that change whole sets of variables, such as the following statement

```python
>>> rc(schema="evolaemp")
```

which switches the variables "asjp", "dolgo", etc. to the ASCII-based transcription system of the ASJP project.
